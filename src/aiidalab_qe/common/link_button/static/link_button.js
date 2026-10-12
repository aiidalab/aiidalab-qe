export const render = ({ model, el }) => {
  const a = document.createElement("a");

  const icon = model.get("icon") || "";
  if (icon) {
    const i = document.createElement("i");
    i.className = `fa fa-${icon}`;
    i.style.marginRight = "0.35em";
    a.appendChild(i);
  }
  const text = document.createTextNode(model.get("description") || "Open");
  a.appendChild(text);

  a.href = model.get("link") || "#";
  a.target = model.get("target") || "_blank";
  const set_tooltip = () => {
    const tooltip = model.get("tooltip") || "";
    a.title = tooltip;
    el.title = tooltip;
  };

  set_tooltip();
  model.on("change:tooltip", set_tooltip);

  const userClasses = (model.get("class_") || "").trim();
  const defaults = "jupyter-button widget-button link-button";
  el.className = [defaults, userClasses].filter(Boolean).join(" ");

  const style_ = model.get("style_") || "";
  if (style_) a.style.cssText = style_;
  else a.removeAttribute("style");

  const set_disabled = () => {
    const disabled = model.get("disabled");
    if (disabled) {
      el.setAttribute("aria-disabled", "true");
      a.setAttribute("aria-disabled", "true");
      a.tabIndex = -1;
    } else {
      el.removeAttribute("aria-disabled");
      a.removeAttribute("aria-disabled");
      a.removeAttribute("tabindex");
    }
  };

  set_disabled();
  model.on("change:disabled", set_disabled);

  a.addEventListener("click", (ev) => {
    if (model.get("disabled")) {
      ev.preventDefault();
      return;
    }

    model.set("clicks", (model.get("clicks") || 0) + 1);
    model.save_changes();

    const prevent = model.get("prevent_default") || !model.get("link");
    if (prevent) {
      ev.preventDefault();
    }
  });

  el.replaceChildren(a);
};
