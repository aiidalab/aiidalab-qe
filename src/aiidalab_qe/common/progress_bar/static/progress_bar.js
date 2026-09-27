export const render = ({ model, el }) => {
  const label = document.createElement("div");
  label.className = "qe-progress-label";

  const track = document.createElement("div");
  track.className = "qe-progress-track";
  track.setAttribute("role", "progressbar");
  track.setAttribute("aria-valuemin", "0");
  track.setAttribute("aria-valuemax", "100");

  const fill = document.createElement("div");
  fill.className = "qe-progress-fill";
  track.appendChild(fill);
  el.classList.add("qe-progress-bar");
  el.replaceChildren(label, track);

  const update_description = () => {
    label.textContent = model.get("description");

    const layout = model.get("description_layout") || {};
    label.style.cssText = "";
    for (const [property, value] of Object.entries(layout)) {
      if (value !== null && value !== "") {
        label.style.setProperty(property.replaceAll("_", "-"), value);
      }
    }
  };

  const update_progress = () => {
    const value = model.get("value");
    const animating = model.get("animating");
    fill.classList.toggle("is-indeterminate", animating);

    if (animating) {
      fill.style.removeProperty("width");
      track.removeAttribute("aria-valuenow");
      track.setAttribute("aria-valuetext", "In progress");
    } else {
      track.setAttribute("aria-valuenow", String(Math.round(value * 100)));
      track.removeAttribute("aria-valuetext");
      fill.style.width = `${value * 100}%`;
    }
  };

  const update_style = () => {
    for (const style of ["info", "warning", "success", "danger", "primary"]) {
      el.classList.toggle(`bar-${style}`, model.get("bar_style") === style);
    }
  };

  update_description();
  update_progress();
  update_style();

  model.on("change:description", update_description);
  model.on("change:description_layout", update_description);
  model.on("change:value", update_progress);
  model.on("change:animating", update_progress);
  model.on("change:bar_style", update_style);

  return () => {
    model.off("change:description", update_description);
    model.off("change:description_layout", update_description);
    model.off("change:value", update_progress);
    model.off("change:animating", update_progress);
    model.off("change:bar_style", update_style);
  };
};
