export const render = ({ model, el }) => {
  let domElement = document.createElement("div");
  el.classList.add("custom-table");
  let selectedIndices = [];

  function drawTable() {
    const data = model.get("data");
    domElement.innerHTML = "";
    let innerHTML =
      "<table><tr>" +
      data[0].map((header) => `<th>${header}</th>`).join("") +
      "</tr>";

    for (let i = 1; i < data.length; i++) {
      innerHTML +=
        "<tr>" + data[i].map((cell) => `<td>${cell}</td>`).join("") + "</tr>";
    }

    innerHTML += "</table>";
    domElement.innerHTML = innerHTML;

    const rows = domElement.querySelectorAll("tr");
    rows.forEach((row, index) => {
      if (index > 0) {
        row.addEventListener("click", () => {
          const rowIndex = index - 1;
          if (selectedIndices.includes(rowIndex)) {
            selectedIndices = selectedIndices.filter((i) => i !== rowIndex);
            row.classList.remove("selected-row");
          } else {
            selectedIndices.push(rowIndex);
            row.classList.add("selected-row");
          }
          model.set("selected_rows", [...selectedIndices]);
          model.save_changes();
        });

        row.addEventListener("mouseover", () => {
          if (!row.classList.contains("selected-row")) {
            row.classList.add("hover-row");
          }
        });

        row.addEventListener("mouseout", () => {
          row.classList.remove("hover-row");
        });
      }
    });
  }

  function updateSelection() {
    const newSelection = model.get("selected_rows");
    selectedIndices = [...newSelection];
    const rows = domElement.querySelectorAll("tr");
    rows.forEach((row, index) => {
      if (index > 0) {
        if (selectedIndices.includes(index - 1)) {
          row.classList.add("selected-row");
        } else {
          row.classList.remove("selected-row");
        }
      }
    });
  }

  drawTable();
  model.on("change:data", drawTable);
  model.on("change:selected_rows", updateSelection);
  el.appendChild(domElement);
};
