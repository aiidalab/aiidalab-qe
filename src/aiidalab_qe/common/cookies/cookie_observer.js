require.undef("cookie_observer");

define("cookie_observer", ["@jupyter-widgets/base"], function (widgets) {
  var CookieObserverView = widgets.DOMWidgetView.extend({
    render: function () {
      const cookieName = this.model.get("name");
      const cookies = document.cookie.split("; ");
      const cookie = cookies.find((row) => row.startsWith(cookieName));
      this.model.set("value", cookie ? true : false);
      this.model.save_changes();
    },
  });

  return {
    CookieObserverView: CookieObserverView,
  };
});
