

from datetime import datetime, timezone

date_time_string = "2026-04-23 16:58:53,814"

dt = datetime.strptime(date_time_string, "%Y-%m-%d %H:%M:%S,%f")
dt = dt.replace(tzinfo=timezone.utc)

unix_ts = dt.timestamp()  # keeps fractional seconds
print(unix_ts)
print(date_time_string)