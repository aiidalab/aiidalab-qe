from __future__ import annotations

import typing as t

import traitlets
from IPython.display import Javascript, display

from .cookie_observer import Cookie


class CookieManager(traitlets.HasTraits):
    """The `CookieManager` class is used to manage browser cookies.

    Attributes
    ----------
    `cookie_name` : `str`
        The name of the cookie.
    `cookie` : `Cookie`
        The cookie observer.
    `cookie_checked` : `traitlets.Bool`
        Whether the cookie has been checked.
    """

    cookies: dict[str, Cookie] = {}

    def __init__(self, cookies: list[Cookie] | None = None, *args, **kwargs):
        """`CookieManager` constructor.

        Parameters
        ----------
        `cookies` : `list[Cookie]`, optional
            One or more cookies to register.
        """
        super().__init__(*args, **kwargs)
        if cookies is not None:
            for cookie in cookies:
                self.add_cookie(cookie)

    def add_cookie(self, cookie_name: str, callback: t.Callable | None = None):
        """Register an observer of a cookie with the given name.

        Parameters
        ----------
        `cookie_name` : `str`
            The name of the cookie.
        `callback` : `callable`, optional
            The callback function to call when the cookie state changes.
        """
        if cookie_name in self.cookies:
            return
        cookie = Cookie(cookie_name)
        cookie.observe(callback, "value")
        self.cookies[cookie.name] = cookie
        display(cookie)

    def has_cookie(self, cookie_name: str) -> bool:
        """Check if the cookie exists.

        Parameters
        ----------
        `cookie_name` : `str`
            The name of the cookie.
        """
        cookie = self.cookies.get(cookie_name, None)
        return cookie is not None and cookie.exists

    def set_cookie(self, cookie_name: str):
        """Set the cookie.

        Parameters
        ----------
        `cookie_name` : `str`
            The name of the cookie.
        """
        js_code = f"""
            (() => {{
                document.cookie = "{cookie_name}=true; path=/";
            }})();
        """
        display(Javascript(js_code))

    def remove_cookie(self, cookie_name: str):
        """Remove the cookie by setting its expiration in the past.

        Parameters
        ----------
        `cookie_name` : `str`
            The name of the cookie.
        """
        past = "Thu, 01 Jan 1970 00:00:00 UTC"
        js_code = f"""
            (() => {{
                document.cookie = "{cookie_name}=; expires={past}; path=/";
            }})();
        """
        display(Javascript(js_code))
