# Security Policy

## Reporting a Vulnerability

If you discover a security vulnerability in jivanu, please report it
responsibly through **GitHub Security Advisories**:

1. Go to the [Security tab](../../security/advisories) of this repository.
2. Click **"Report a vulnerability"**.
3. Fill in the details and submit.

**Do not open a public issue for security vulnerabilities.**

## Response Timeline

| Action | Target |
|---|---|
| Acknowledgement | Within **48 hours** |
| Initial assessment | Within **5 business days** |
| Fix for critical severity | Within **14 days** |
| Fix for high severity | Within **30 days** |
| Fix for moderate/low severity | Next scheduled release |

## Supported Versions

| Version | Supported |
|---------|-----------|
| 1.x | Yes |
| < 1.0 | No |

## Scope

This policy covers the `jivanu` crate and its published API. Vulnerabilities
in dependencies (serde, thiserror, hisab, etc.) should be reported to the
respective upstream projects and flagged here if they affect jivanu users.

## Design Principles

- **No unwrap/panic** — all error paths return `Result` or safe defaults
- **No unsafe** — zero unsafe blocks in library code
- **No I/O** — pure computation, no network or filesystem access
- **Validated inputs** — all public functions validate parameters before use
- **Minimal dependencies** — serde, thiserror, tracing (required); hisab, kimiya, pramana (optional)
- **Supply chain** — `cargo audit` and `cargo deny` in CI on every push

## Disclosure

We follow coordinated disclosure. Once a fix is released, we will publish a
security advisory crediting the reporter (unless anonymity is requested).
