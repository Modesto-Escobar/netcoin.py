"""Timelines of periods and events (ports of ``timeline_rd3``, ``timeCoin``, ``incTime``, ``mobileEdges``)."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from .similarity import edge_list
from .webgraph import WebGraph, check_cex, check_language, records


def _as_seconds(series: pd.Series) -> pd.Series:
    """Date-times as seconds since 1970-01-01 (R's ``as.numeric`` of a POSIXct)."""
    stamps = pd.to_datetime(series)
    if getattr(stamps.dt, "tz", None) is not None:
        stamps = stamps.dt.tz_convert("UTC").dt.tz_localize(None)
    return (stamps - pd.Timestamp("1970-01-01")).dt.total_seconds()


class TimeCoin(WebGraph):
    """A timeline: bars for periods (people, schools...), optionally with events on them.

    Port of rD3plot's ``timeline_rd3``. Build it with :func:`time_coin`.

    Attributes
    ----------
    periods : DataFrame
        One row per period, with its name, start and end.
    events : DataFrame or None
        Events placed on the periods.
    options : dict
        What the page reads.
    """

    _engine = "timeline.js"

    def __init__(
        self,
        periods: pd.DataFrame,
        name: str = "name",
        start: str = "start",
        end: str = "end",
        group: str | None = None,
        text: str | None = None,
        main: str | None = None,
        note: str | None = None,
        info: str | None = None,
        events: pd.DataFrame | None = None,
        event_names: str = "name",
        event_period: str = "period",
        event_time: str = "date",
        event_color: str | None = None,
        event_shape: str | None = None,
        cex: float = 1,
        language: str = "en",
        directory: str | None = None,
    ):
        missing = {name, start, end} - set(periods.columns)
        if missing:
            raise ValueError("name, start and end: must be present in periods data frame as columns.")
        periods = periods.copy()

        options = {"name": name, "start": start, "end": end, "cex": check_cex(cex)}
        for key, value in (("group", group), ("text", text), ("main", main), ("note", note),
                           ("info", info)):
            if value is not None:
                options[key] = value
        options["language"] = check_language(language)

        if pd.api.types.is_datetime64_any_dtype(periods[start]) and \
                pd.api.types.is_datetime64_any_dtype(periods[end]):
            options["POSIXct"] = True
        for column in (start, end):
            if pd.api.types.is_datetime64_any_dtype(periods[column]):
                periods[column] = _as_seconds(periods[column])

        if events is not None:
            events = events.copy()
            if pd.api.types.is_datetime64_any_dtype(events[event_time]):
                events[event_time] = _as_seconds(events[event_time])
            events = events.sort_values(event_time, kind="stable").reset_index(drop=True)
            options["eventNames"] = event_names
            options["eventPeriod"] = event_period
            options["eventTime"] = event_time
            if event_color is not None:
                options["eventColor"] = event_color
            if event_shape is not None:
                options["eventShape"] = event_shape

        self.periods = periods
        self.events = events
        self.options = options
        if directory:
            self._create(directory)

    def _data(self) -> dict:
        data = {"periods": records(self.periods)}
        if self.events is not None:
            data["events"] = records(self.events)
        data["options"] = self.options
        return data

    def summary(self):
        """Print how many periods there are and the distribution of start and end."""
        print(f"{len(self.periods)} categories.")
        for key in ("start", "end"):
            column = self.options[key]
            print(f"{column}'s distribution:")
            print(self.periods[column].describe()[["min", "25%", "50%", "mean", "75%", "max"]]
                  .round(2).to_string())


def time_coin(
    periods: pd.DataFrame,
    name: str = "name",
    start: str = "start",
    end: str = "end",
    group: str | None = None,
    text: str | None = None,
    main: str | None = None,
    note: str | None = None,
    info: str | None = None,
    events: pd.DataFrame | None = None,
    event_names: str = "name",
    event_period: str = "period",
    event_time: str = "date",
    event_color: str | None = None,
    event_shape: str | None = None,
    cex: float = 1,
    language: str = "en",
    directory: str | None = None,
) -> TimeCoin:
    """Timeline of periods (and events) from a table.

    Port of ``timeCoin`` from the R package netCoin.

    Parameters
    ----------
    periods : DataFrame
        One row per period: ``name``, ``start`` and ``end`` (numbers or date-times).
    group, text, info
        Columns used to colour the bars, to describe them and for the info panel.
    events : DataFrame, optional
        Events: the ``event_names`` column names them, ``event_period`` says which
        period they belong to (a name of ``periods``) and ``event_time`` when
        (a number or a date-time).
    event_color, event_shape
        Columns of ``events`` that set the colour and the shape of the marks.
    main, note, cex, language, directory
        Title, footnote, font size, ``"en"``/``"es"``/``"ca"`` and where to write the page.

    Returns
    -------
    TimeCoin
    """
    return TimeCoin(periods, name, start, end, group, text, main, note, info, events, event_names,
                    event_period, event_time, event_color, event_shape, cex, language, directory)


def inc_time(data: pd.DataFrame, name: str = "name", beginning: str = "birth",
             end: str = "death") -> pd.DataFrame:
    """Which of the people (or periods) are alive in every year.

    Port of ``incTime``: a table with one row per year, from the first beginning to the
    last end, and one 0/1 column per person. An open end (missing) goes on forever.
    A missing beginning gives missing values.
    """
    start = pd.to_numeric(data[beginning], errors="coerce").to_numpy(dtype=float)
    stop = pd.to_numeric(data[end], errors="coerce").to_numpy(dtype=float)
    years = np.arange(int(np.nanmin(start)), int(np.nanmax(stop)) + 1)
    inside = (years[:, None] >= start[None, :]) & ((years[:, None] <= stop[None, :]) | np.isnan(stop)[None, :])
    values = inside.astype(float)
    values[:, np.isnan(start)] = np.nan
    out = pd.DataFrame(values, index=years, columns=[str(v) for v in data[name]])
    return out.astype("Int64") if np.isnan(values).any() else out.astype(int)


def mobile_edges(data: pd.DataFrame, name: int | str = 0, number: int | str = 1,
                 difference: float = 0) -> pd.DataFrame:
    """Links between the items whose numbers are close (``|a - b| <= difference``).

    Port of ``mobileEdges``: ``sim.`` is ``1 + difference - |a - b|`` (at least 1) and
    ``dist.`` its complement ``1 + difference - sim.``. Use it, for instance, to link
    people born in the same period.
    """
    labels = data.iloc[:, name] if isinstance(name, int) else data[name]
    values = pd.to_numeric(data.iloc[:, number] if isinstance(number, int) else data[number],
                           errors="coerce").to_numpy(dtype=float)
    gap = np.abs(values[:, None] - values[None, :])
    similarity = np.where(gap <= difference, 1 + difference - gap, 0.0)
    np.fill_diagonal(similarity, 0.0)
    matrix = pd.DataFrame(similarity, index=[str(v) for v in labels], columns=[str(v) for v in labels])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        links = edge_list(matrix, "shape", min=1)
    links = links.rename(columns={"value": "sim."})
    links["dist."] = 1 + difference - links["sim."]
    return links.reset_index(drop=True)
