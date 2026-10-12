"""Model client for Alibaba Cloud Model Studio, talking the OpenAI-compatible
chat-completions protocol over the standard library only.

Credential handling
-------------------
The key is read from an environment variable at call time. It is never written
to a config file, a source file, a run record, a log or a report. Every string
that could be logged passes through ``redact()`` first, so an accidental echo
of an Authorization header cannot leak the value. If the variable is missing the
client refuses to construct -- there is no fallback, no prompt, no default.

Transport policy
----------------
Only transport failures are retried: timeouts, connection errors, HTTP 429 and
5xx. A response that came back but contains a wrong answer is NEVER retried --
re-rolling until a good answer appears would invalidate the comparison. The
retry count and backoff come from run_config.json and are recorded per request.
"""
from __future__ import annotations

import json
import os
import re
import time
import urllib.error
import urllib.request

RETRYABLE_STATUS = {408, 409, 425, 429, 500, 502, 503, 504}

# Credential variables for THIS EXPERIMENT, in preference order. They are
# deliberately separate from anything the Claude Code harness uses, so the
# experiment can never silently spend the coding agent's subscription, and so
# revoking one does not break the other.
EXPERIMENT_API_KEY_ENVS = ('BMB_API_KEY', 'DASHSCOPE_API_KEY')

# --------------------------------------------------------------------------
# providers
# --------------------------------------------------------------------------
#
# The call channel is an explicit choice, not a hard-coded constant. Two are
# defined; the default stays `bailian`, so existing configurations and the
# standard path are unaffected by anything below.
#
# The Token Plan entry carries a verbatim quotation of the platform's stated
# scope. Selecting it does NOT make the use permitted, and nothing in this
# project claims otherwise: the warning is printed at startup and written into
# the run manifest so the choice is visible in every artifact.

TOKEN_PLAN_TERMS_CN = (
    '仅限在编程工具和智能体工具（如 Claude Code、Cursor、Qwen Code、Qoder、Qoder CN、'
    'OpenClaw 等）中交互式使用，不可用于自动化脚本、自定义应用程序后端或任何非交互式'
    '批量调用场景。……将套餐 API Key 用于允许范围之外的调用将被视为违规或滥用，'
    '可能会导致订阅被暂停或 API Key 被封禁。')

TOKEN_PLAN_TERMS_SOURCE = (
    'https://help.aliyun.com/zh/model-studio/token-plan-personal-overview'
    '（订阅前须知 / 使用范围，本机 2026-10-06 抓取，HTTP 200）')

PROVIDERS = {
    'bailian': {
        'label': 'Model Studio standard model API (pay-as-you-go)',
        'label_cn': '百炼普通模型 API（按量计费）',
        'credential_envs': ('BMB_API_KEY', 'DASHSCOPE_API_KEY'),
        'default_base_url': 'https://dashscope.aliyuncs.com/compatible-mode/v1',
        'workspace_scoped_hosts_allowed': False,
        'terms_warning': None,
    },
    'token_plan': {
        'label': 'Token Plan (personal subscription)',
        'label_cn': 'Token Plan（个人版，订阅套餐）',
        'credential_envs': ('BMB_TOKEN_PLAN_API_KEY',),
        'default_base_url': 'https://token-plan.cn-beijing.maas.aliyuncs.com/compatible-mode/v1',
        'workspace_scoped_hosts_allowed': True,
        'terms_warning': (
            'The platform states this plan is for interactive use inside coding and agent '
            'tools and forbids "automation scripts, custom application backends, or any '
            'non-interactive batch calling scenario", warning that out-of-scope use may '
            'suspend the subscription or ban the key. THIS PROJECT DOES NOT CLAIM THIS USE '
            'IS PERMITTED. The user selected this provider explicitly and knowingly; the '
            'choice and this warning are recorded in every run manifest.'),
        'terms_verbatim_cn': TOKEN_PLAN_TERMS_CN,
        'terms_source': TOKEN_PLAN_TERMS_SOURCE,
    },
}

# Hosts that are never valid for the standard provider. Kept from the earlier
# implementation so the standard path cannot silently drift onto a plan or
# workspace-scoped host.
BLOCKED_FOR_STANDARD_PROVIDER = {
    'token-plan.': ('a Token Plan host. Set provider=token_plan to use it deliberately, with '
                    'the matching credential variable.'),
    '.maas.aliyuncs.com': ('a workspace-scoped MaaS host, which is tied to a specific '
                           'workspace subscription; the standard provider uses the DashScope '
                           'OpenAI-compatible host.'),
}


class ModelConfigMissing(RuntimeError):
    """Raised before a run starts when the endpoint is not configured."""


class EndpointNotPermitted(RuntimeError):
    """Raised when the configured endpoint is not licensed for this use."""


def check_endpoint_permitted(base_url, provider='bailian'):
    """Validate the endpoint against the SELECTED provider.

    This is a consistency check, not a claim of permission. For the standard
    provider it refuses plan and workspace-scoped hosts outright. For the
    Token Plan provider it allows that host but never asserts the use is
    sanctioned -- the terms warning travels with the provider and is recorded
    in the run manifest.
    """
    if not base_url:
        return
    if provider not in PROVIDERS:
        raise ModelConfigMissing(
            'Unknown provider %r. Known providers: %s'
            % (provider, ', '.join(sorted(PROVIDERS))))
    lowered = base_url.lower()
    if provider == 'bailian':
        for marker, reason in BLOCKED_FOR_STANDARD_PROVIDER.items():
            if marker in lowered:
                raise EndpointNotPermitted(
                    'Endpoint %s does not belong to provider "bailian".\n%s\n'
                    'The standard host is https://dashscope.aliyuncs.com/compatible-mode/v1'
                    % (base_url, reason))
        if not lowered.startswith('https://'):
            raise EndpointNotPermitted('Endpoint must use HTTPS: %s' % base_url)
        return
    # token_plan: the plan host is expected; a mismatch is almost certainly a
    # configuration slip and would fail authentication anyway
    if 'token-plan.' not in lowered and '.maas.aliyuncs.com' not in lowered:
        raise EndpointNotPermitted(
            'provider=token_plan but the endpoint does not look like a Token Plan host: %s\n'
            'Expected something like %s'
            % (base_url, PROVIDERS['token_plan']['default_base_url']))
    if not lowered.startswith('https://'):
        raise EndpointNotPermitted('Endpoint must use HTTPS: %s' % base_url)


def provider_terms_warning(provider):
    """The warning to print and to persist for the selected provider."""
    spec = PROVIDERS.get(provider) or {}
    return spec.get('terms_warning')


def credential_env_names():
    """Every environment variable that may hold a secret for this project."""
    names = list(EXPERIMENT_API_KEY_ENVS)
    for name in os.environ:
        upper = name.upper()
        if upper in names:
            continue
        if upper.endswith('API_KEY') or upper.endswith('_TOKEN') or 'SECRET' in upper:
            names.append(name)
    return names


def resolve_api_key(explicit_env_names=None):
    """The experiment key, taken from the environment. Never logged."""
    for name in (explicit_env_names or EXPERIMENT_API_KEY_ENVS):
        value = os.environ.get(name)
        if value:
            return value, name
    return None, None


def redact(text) -> str:
    """Mask anything that looks like a bearer token or an API key.

    Implemented as a single forward pass. An earlier version replaced in a loop
    and matched its own output forever; the selftest hangs if that returns.
    """
    if text is None:
        return ''
    text = str(text)
    for name in credential_env_names():
        value = os.environ.get(name)
        if value and len(value) >= 8:
            text = text.replace(value, '<redacted>')
    text = re.sub(r'(?i)\bBearer\s+[^\s"\',;)]+', 'Bearer <redacted>', text)
    text = re.sub(r'\bsk-[A-Za-z0-9_\-]{6,}', 'sk-<redacted>', text)
    return text


class ModelClient:
    """One endpoint. Immutable for the life of a run."""

    def __init__(self, base_url, model_id, api_key_env='DASHSCOPE_API_KEY',
                 temperature=0.0, top_p=1.0, max_tokens=512,
                 timeout=120, max_attempts=3, backoff_seconds=2.0,
                 use_json_object=True, region=None, extra_body=None,
                 provider='bailian'):
        if not base_url or not model_id:
            raise ModelConfigMissing(
                'base_url and model_id must both be set. Configure them in '
                'config/model_config.json under "active", or export '
                'BMB_BASE_URL / BMB_MODEL_ID. See README.md.')
        self.provider = provider
        check_endpoint_permitted(base_url, provider)
        self.api_key_env = api_key_env
        self.base_url = base_url.rstrip('/')
        self.model_id = model_id
        self.region = region
        self.temperature = temperature
        self.top_p = top_p
        self.max_tokens = max_tokens
        self.timeout = timeout
        self.max_attempts = max_attempts
        self.backoff_seconds = backoff_seconds
        self.use_json_object = use_json_object
        # Provider-specific request fields, e.g. {"enable_thinking": false}.
        # These are NOT part of the OpenAI schema, so they are kept separate
        # and recorded in the run manifest rather than hard-coded here.
        self.extra_body = dict(extra_body or {})

    # -- credentials -------------------------------------------------------
    def credential_envs(self):
        """Credential variables for the SELECTED provider, most specific first.

        token_plan looks only at BMB_TOKEN_PLAN_API_KEY: the plan credential is
        never borrowed from a shared variable, and nothing here reads another
        tool's credential file.
        """
        provider_envs = tuple(PROVIDERS.get(self.provider, {}).get('credential_envs') or ())
        if self.api_key_env:
            provider_envs = (self.api_key_env,) + provider_envs
        return tuple(dict.fromkeys(provider_envs))

    @property
    def api_key(self) -> str:
        key, _name = resolve_api_key(self.credential_envs())
        if not key:
            raise ModelConfigMissing(
                'No API key found for provider %r. Set %s in this shell only. '
                'Do not put it in a file this project tracks, and do not paste it into '
                'a chat. See README.md for the local configuration steps.'
                % (self.provider, ' or '.join(self.credential_envs())))
        return key

    @property
    def api_key_var_in_use(self):
        return resolve_api_key(self.credential_envs())[1]

    def credential_present(self) -> bool:
        return resolve_api_key(self.credential_envs())[0] is not None

    # -- request -----------------------------------------------------------
    def _endpoint(self):
        return self.base_url + '/chat/completions'

    def _payload(self, messages):
        body = {
            'model': self.model_id,
            'messages': messages,
            'temperature': self.temperature,
            'top_p': self.top_p,
            'max_tokens': self.max_tokens,
            'stream': False,
        }
        if self.use_json_object:
            body['response_format'] = {'type': 'json_object'}
        body.update(self.extra_body)
        return body

    def chat(self, messages):
        """Return a dict with the response text, usage, latency and attempt log."""
        payload = json.dumps(self._payload(messages), ensure_ascii=False).encode('utf-8')
        headers = {
            'Content-Type': 'application/json',
            'Authorization': 'Bearer ' + self.api_key,
        }
        attempts = []
        last_error = 'unknown'
        for attempt in range(1, self.max_attempts + 1):
            started = time.time()
            try:
                request = urllib.request.Request(self._endpoint(), data=payload,
                                                 headers=headers, method='POST')
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    raw = response.read().decode('utf-8', errors='replace')
                latency = time.time() - started
                try:
                    parsed = json.loads(raw)
                except json.JSONDecodeError as exc:
                    return {'ok': False, 'transport_error': 'malformed response body: %s' % exc,
                            'attempts': attempts, 'latency_s': latency}
                text, usage = self._extract(parsed)
                return {'ok': True, 'text': text, 'usage': usage, 'latency_s': latency,
                        'attempts': attempts, 'finish_reason': self._finish_reason(parsed)}
            except urllib.error.HTTPError as exc:
                latency = time.time() - started
                body = ''
                try:
                    body = exc.read().decode('utf-8', errors='replace')
                except Exception:  # noqa: BLE001
                    pass
                detail = 'HTTP %s: %s' % (exc.code, redact(body)[:400])
                attempts.append({'attempt': attempt, 'outcome': detail,
                                 'latency_s': latency, 'retryable': exc.code in RETRYABLE_STATUS})
                last_error = detail
                if exc.code not in RETRYABLE_STATUS:
                    return {'ok': False, 'transport_error': detail, 'attempts': attempts,
                            'latency_s': latency}
            except Exception as exc:  # noqa: BLE001 - timeouts and connection errors
                latency = time.time() - started
                detail = '%s: %s' % (type(exc).__name__, redact(exc))
                attempts.append({'attempt': attempt, 'outcome': detail,
                                 'latency_s': latency, 'retryable': True})
                last_error = detail
            if attempt < self.max_attempts:
                time.sleep(self.backoff_seconds * attempt)
        return {'ok': False, 'transport_error': last_error, 'attempts': attempts,
                'latency_s': None}

    @staticmethod
    def _extract(parsed):
        try:
            choice = parsed['choices'][0]
        except (KeyError, IndexError, TypeError):
            return '', None
        message = choice.get('message') or {}
        text = message.get('content')
        if text is None:
            text = choice.get('text') or ''
        return text, parsed.get('usage')

    @staticmethod
    def _finish_reason(parsed):
        try:
            return parsed['choices'][0].get('finish_reason')
        except (KeyError, IndexError, TypeError):
            return None


def client_from_config(model_config, run_config, model_id=None, base_url=None,
                       region=None):
    """Build a client from the model config, honouring environment overrides.

    Nothing is guessed: if the active profile is unset and the environment does
    not override it, ModelConfigMissing is raised and the run does not start.
    """
    active = dict(model_config.get('active') or {})
    env = model_config.get('env_overrides') or {}
    provider = (active.get('provider') or 'bailian')
    region = region or os.environ.get(env.get('region', ''), '') or active.get('region')
    base_url = base_url or os.environ.get(env.get('base_url', ''), '') or active.get('base_url')
    model_id = model_id or os.environ.get(env.get('model_id', ''), '') or active.get('model_id')

    if region and not base_url:
        region_block = (model_config.get('regions') or {}).get(region) or {}
        base_url = region_block.get('base_url_legacy') or region_block.get('base_url_workspace_scoped')
        if base_url and '{WorkspaceId}' in base_url:
            base_url = None  # cannot be resolved without the user's workspace id

    generation = run_config.get('generation', {})
    retry = run_config.get('retry', {})
    return ModelClient(
        base_url=base_url, model_id=model_id,
        api_key_env=model_config.get('api_key_env', 'DASHSCOPE_API_KEY'),
        temperature=generation.get('temperature', 0.0),
        top_p=generation.get('top_p', 1.0),
        max_tokens=generation.get('max_tokens', 512),
        max_attempts=retry.get('max_attempts', 3),
        backoff_seconds=retry.get('backoff_seconds', 2.0),
        region=region,
        extra_body=generation.get('extra_body') or {},
        provider=provider,
    )


def credential_env_names_for(provider):
    return tuple(PROVIDERS.get(provider, {}).get('credential_envs') or ())