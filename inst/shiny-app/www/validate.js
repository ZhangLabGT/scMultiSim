function setInvalid(el, msg) {
  el.classList.add("is-invalid");
  let sb = el.nextElementSibling;
  while (sb !== null && !sb.classList.contains("invalid-feedback")) {
    sb = sb.nextElementSibling;
  }
  if (sb === null) {
    console.error("No sibling with class 'invalid-feedback' found");
  }
  sb.innerText = msg;
}

function setValid(el) {
  el.classList.remove("is-invalid");
}

function validateMultipleElements(els, custom_fn) {
  const [valid, msg] = custom_fn(els);
  if (valid) {
    els.forEach((el) => setValid(el));
  } else {
    els.forEach((el) => setInvalid(el, msg));
  }
}

function validateSingleElement(el, custom_fn) {
  if (typeof custom_fn === "function") {
    const [valid, msg] = custom_fn(el);
    if (valid) {
      setValid(el);
    } else {
      setInvalid(el, msg);
    }
    return;
  }

  const data = el.dataset.v;
  const [type, valueRange] = data.split(":");

  let value = el.value;
  if (type === "i" || type === "f") {
    if (type === "i" && /[^\d]/.test(value)) {
      setInvalid(el, "Please enter an integer number");
      return;
    }
    value = type === "i" ? parseInt(value) : parseFloat(value);
    if (isNaN(value)) {
      setInvalid(el, `Please enter a number`);
      return;
    }
  } else if (type === "n") {
    return;
  }

  if (valueRange !== undefined) {
    const [min, max] = valueRange.split("-");
    if (value < min || value > max) {
      setInvalid(el, `Please enter a value between ${min} and ${max}`);
      return;
    }
  }

  setValid(el);
}

function validate(name, custom_fn) {
  if (Array.isArray(name)) {
    const els = name.map((n) => document.querySelector(`[name="${n}"]`));
    validateMultipleElements(els, custom_fn);
    for (const el of els) {
      el.addEventListener("input", () => {
        validateMultipleElements(els, custom_fn);
      });
    }
    return;
  }
  const el = document.querySelectorAll(`[name="${name}"]`);
  if (el.length === 0) {
    console.error(`Element with name ${name} not found`);
  } else if (el.length === 1) {
    validateSingleElement(el[0], custom_fn);
    el[0].addEventListener("input", () => {
      validateSingleElement(el[0], custom_fn);
    });
  }
}

function init() {
  const validatedNames = new Set();

  function customValidate(n, fn) {
    validate(n, fn);
    validatedNames.add(n);
  }

  customValidate("discretePopSize", (el) => {
    const popType = document.querySelector("input[name='_cellPop']:checked").value;
    if (popType === "continuous" || el.value.length === 0) return [true, null];
    const nClus = parseInt(document.querySelector("input[name='_numCluster']").value);
    const nCell = parseInt(document.querySelector("input[name='num.cells']").value);
    const values = el.value.split(",").map((v) => parseInt(v));
    if (values.length !== nClus) return [false, "Number of clusters not matching"];
    if (values.reduce((a, b) => a + b, 0) !== nCell) return [false, "Number of cells not matching"];
    return [true, null];
  });

  customValidate(["_regionDist0", "_regionDist1", "_regionDist2"], (el) => {
    const distValues = [0, 1, 2].map((i) => parseFloat(document.querySelector(`input[name='_regionDist${i}']`).value));
    const currValue = parseFloat(el.value);
    if (currValue < 0 || currValue > 1) return [false, "Please enter a value between 0 and 1"];
    if (distValues.reduce((a, b) => a + b, 0) !== 1) return [false, "Sum of all values should be 1"];
    return [true, null];
  });

  const allInputs = document.querySelectorAll("input");
  for (const input of allInputs) {
    if (typeof input.dataset.v === "string") {
      if (input.dataset.v !== "c" && !validatedNames.has(input.name)) {
        validate(input.name);
        validatedNames.add(input.name);
      }
    } else {
      if (!(input.type === "radio" || input.type === "checkbox")) {
        console.warn("Validation: data attribute not found", input);
      }
    }
  }
}

document.addEventListener("DOMContentLoaded", init);
