from datetime import datetime

from dateutil.relativedelta import relativedelta


def format_time(time: datetime):
    return time.strftime("%Y-%m-%d %H:%M:%S")


def relative_time(time: datetime):
    # TODO consider using humanize or arrow libraries for this
    now = datetime.now(time.tzinfo)
    delta = relativedelta(now, time)
    if delta.years > 0:
        return f"{delta.years} year{'s' if delta.years > 1 else ''} ago"
    elif delta.months > 0:
        return f"{delta.months} month{'s' if delta.months > 1 else ''} ago"
    elif delta.days > 0:
        return f"{delta.days} day{'s' if delta.days > 1 else ''} ago"
    elif delta.hours > 0:
        return f"{delta.hours} hour{'s' if delta.hours > 1 else ''} ago"
    elif delta.minutes > 0:
        return f"{delta.minutes} minute{'s' if delta.minutes > 1 else ''} ago"
    else:
        return f"{delta.seconds} second{'s' if delta.seconds > 1 else ''} ago"
