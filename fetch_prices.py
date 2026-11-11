"""Download and freeze Yahoo Finance daily prices used for Assignment 2."""

# The report runs from the CSVs already in raw/. This file records how they
# were obtained; rerunning it may replace those CSVs with revised Yahoo data.

import csv
import datetime as dt
import json
import urllib.parse
import urllib.request
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT = Path(__file__).parent
RAW = ROOT / "raw"
RAW.mkdir(exist_ok=True)
SYMBOLS = ["INTC", "SSAB-A.ST", "SEK=X"]
START = dt.datetime(2009, 10, 1, tzinfo=dt.timezone.utc)
# The API's period2 is exclusive. The last required share close is 30 June.
END = dt.datetime(2026, 7, 1, tzinfo=dt.timezone.utc)

for symbol in SYMBOLS:
    query = urllib.parse.urlencode(
        {
            "period1": int(START.timestamp()),
            "period2": int(END.timestamp()),
            "interval": "1d",
            "events": "div,splits",
            "includeAdjustedClose": "true",
        }
    )
    url = f"https://query1.finance.yahoo.com/v8/finance/chart/{urllib.parse.quote(symbol)}?{query}"
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = json.load(response)
    if payload["chart"]["error"] is not None:
        raise RuntimeError((symbol, payload["chart"]["error"]))
    result = payload["chart"]["result"][0]
    with (RAW / f"{symbol.replace('=', '_')}.json").open("w") as stream:
        json.dump(payload, stream)
    # A daily timestamp is the exchange's trading date, which can differ
    # from the UTC date, especially for the USD/SEK series.
    timezone = ZoneInfo(result["meta"]["exchangeTimezoneName"])
    quote = result["indicators"]["quote"][0]
    adjusted = result["indicators"].get("adjclose", [{}])[0].get("adjclose", [])
    # Keep Close for the assignment's price-only model and Adj Close for the
    # separate distribution-adjusted comparison in the report.
    records = []
    for i, stamp in enumerate(result["timestamp"]):
        records.append(
            {
                "date": dt.datetime.fromtimestamp(stamp, timezone).date().isoformat(),
                "close": quote["close"][i],
                "adjusted_close": adjusted[i] if adjusted else None,
                "currency": result["meta"]["currency"],
            }
        )
    with (RAW / f"{symbol.replace('=', '_')}.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(symbol, len(records), records[0], records[-1], "events", list(result.get("events", {})))

# Record the retrieval time so the frozen files can be traced later.
with (RAW / "provenance.json").open("w") as stream:
    json.dump(
        {
            "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "source": "Yahoo Finance chart API; Yahoo Finance historical pages listed in report",
            "symbols": SYMBOLS,
            "requested_start": START.date().isoformat(),
            "requested_end_exclusive": END.date().isoformat(),
            "field_used": "close",
        },
        stream,
        indent=2,
    )
