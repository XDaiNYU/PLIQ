# GitHub upload (clean package)

**Repository:** https://github.com/XDaiNYU/PLIQ

## What is included

- Source: `pliq/`, `pliq_edockq_scoring.py`, `examples/`, `scripts/`, `tests/`
- Docs: `README.md`, `docs/`, `reports/`, `CHANGELOG.md`, `LICENSE`
- CI: `.github/workflows/ci.yml`
- Git hook: `scripts/git-hooks/prepare-commit-msg` (strips Cursor co-author from commits)

## Full replace `main` (v0.6.1)

```bash
cd pliq_ver6
git init
cp scripts/git-hooks/prepare-commit-msg .git/hooks/prepare-commit-msg
chmod +x .git/hooks/prepare-commit-msg
git add .
git commit -m "Release PLIQ v0.6.1"
git remote add origin https://github.com/XDaiNYU/PLIQ.git
git branch -M main
git push -u origin main --force
```

## No Cursor in commits

1. Install the hook above before `git commit`.
2. In Cursor: **Settings → Agent → Attribution → disable Commit Attribution**.

## No Cursor as GitHub collaborator

GitHub → **Settings → Collaborators** on `XDaiNYU/PLIQ` → remove any Cursor-related account if present.
