"""
Entry point — run with:
  uvicorn sl_agent.main:app --reload --port 8000
or
  python -m sl_agent.main
"""
import uvicorn

from sl_agent.api.app import app  # noqa: F401 — re-exported for uvicorn

if __name__ == "__main__":
    uvicorn.run(
        "sl_agent.main:app",
        host="0.0.0.0",
        port=8000,
        reload=False,
        workers=1,
        log_level="info",
    )
