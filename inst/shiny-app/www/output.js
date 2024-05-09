const defaultValues = {};

function outputOptions() {
  const options = {};

  function add(name, getter = null, inputName = null) {
    const defaultValue = defaultValues[name];
    let value;
    if (getter) {
      value = typeof getter === "function" ? getter() : getter;
    } else {
      const el = document.querySelector(`[name="${inputName || name}"]`);
      if (el === null) {
        console.error(`Element with name ${name} not found`);
        return;
      }
      value = el.type === "checkbox" ? el.checked.toString().toUpperCase() : el.value;
    }
    if (value !== defaultValue) {
      options[name] = value;
    }
  }

  const isDiscrete = document.querySelector("[name=_cellPop]:checked").value === "discrete";
  const grnEnabled = document.querySelector("[name=_grnEnabled]:checked").value === "on";

  add("seed");
  add("speed.up");

  add("num.cifs");
  add("diff.cif.fraction");
  add("cif.center");
  add("cif.sigma");
  add("giv.mean");
  add("giv.sd");
  add("giv.prob");
  add("num.cells");

  if (isDiscrete) {
    add("discrete.cif", "TRUE");
    const numCluster = document.querySelector("[name=_numCluster]").value;
    add("tree", `rtree(${numCluster})`);
    add("discrete.min.pop.size");
    const discretePopSize = document.querySelector("[name=discretePopSize]").value;
    if (discretePopSize.length > 0) {
      add("discrete.pop.size", `as.integer(c(${discretePopSize}))`);
    }
  } else {
    add("tree", () => {
      let treeVal = document.querySelector("[name=tree]:checked").value;
      if (treeVal === "custom") {
        treeVal = document.querySelector("[name='_treeCustom']").value;
      } else {
        treeVal = {
          phyla1: "Phyla1()",
          phyla3: "Phyla3()",
          phyla5: "Phyla5()",
        }[treeVal];
      }
      return treeVal;
    });
    add("use.impulse");
  }

  if (grnEnabled) {
    add("GRN");
    add("num.genes");
  } else {
    add("GRN", "NA");
    add("num.genes", null, "num.genes2");
  }

  const useCustomATACDensity = document.querySelector("[name=useCustomATACDensity]").checked;

  add("riv.mean");
  add("riv.sd");
  add("riv.prob");

  add("region.distrib", () => {
    const v = [0, 1, 2].map((i) => document.querySelector(`[name=_regionDist${i}]`).value).join(", ");
    return `c(${v})`;
  });
  add("atac.effect");
  add("atac.p_zero");
  if (useCustomATACDensity) {
    add("atac.density");
  }

  add("scale.s");
  add("bimod");

  const velocityEnabled = document.querySelector("[name=velocity]:checked").value === "on";
  if (velocityEnabled) {
    add("do.velocity", "TRUE");
    add("beta");
    add("d");
    add("num.cycles");
    add("cycle.len");
  } else {
    add("intrinsic.noise");
  }

  let cciOptions = null;
  function add_sp(name, getter = null) {
    let value;
    if (getter) {
      value = typeof getter === "function" ? getter() : getter;
    } else {
      const el = document.querySelector(`[name="sp.${name}"]`);
      if (el === null) {
        console.error(`Element with name ${name} not found`);
        return;
      }
      value = el.type === "checkbox" ? el.checked.toString().toUpperCase() : el.value;
    }
    cciOptions[name] = value;
  }

  const spatialEnabled = document.querySelector("[name=_spatialEnabled]").checked;
  if (spatialEnabled) {
    cciOptions = {};

    add_sp("grid.size");
    add_sp("step.size");
    add_sp("max.neighbors");
    add_sp("params", () => document.querySelector("[name=_CCI]").value);

    const spLayout = document.querySelector("[name=_spLayout]").value;
    if (spLayout === "normal") {
      add_sp("layout", `"enhanced"`);
      add_sp("same.type.prob");
    } else if (spLayout === "islands") {
      const spIslands = document.querySelector("[name=_spIslands]").value;
      add_sp("layout", `"islands:${spIslands}"`);
    } else if (spLayout === "layers") {
      add_sp("layout", `"layers"`);
    }

    if (cciOptions.params.length === 0) {
      alert("Please select a CCI dataframe to simulate spatial data.");
    }
  }

  console.log(options, cciOptions);

  const optString = `options <- list(
    ${Object.entries(options)
      .map(([k, v]) => `${k} = ${v}`)
      .join(", ")}
    ${
      cciOptions
        ? `, cci = list(${Object.entries(cciOptions)
            .map(([k, v]) => `${k} = ${v}`)
            .join(", ")})`
        : ""
    }
    )`;
  Shiny.setInputValue("generatedOptions", optString);
  const modal = new bootstrap.Modal(document.getElementById("outputModal"));
  modal.show();
}

function init() {
  const outputButton = document.querySelector("#outputButton");
  outputButton.addEventListener("click", outputOptions);

  document.querySelector("#resetButton").addEventListener("click", () => {
    location.reload();
  });

  document.querySelector("#stopApp").addEventListener("click", () => {
    Shiny.setInputValue("stopApp", "YES");
    close();
  });

  Shiny.addCustomMessageHandler("Defaults", (v) => {
    console.log("Def", v);
    Object.assign(defaultValues, v);
  });
}

document.addEventListener("DOMContentLoaded", init);
