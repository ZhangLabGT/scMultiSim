const colorList = [
  "#e60049",
  "#0bb4ff",
  "#50e991",
  "#e6d800",
  "#9b19f5",
  "#ffa300",
  "#dc0ab4",
  "#b3d4ff",
  "#00bfa0",
  "#ff6f00",
  "#ff00b4",
  "#00ff6f",
  "#b4ff00",
  "#6f00ff",
  "#00ff00",
  "#ff0000",
  "#0000ff",
  "#00ffff",
  "#ffff00",
];

const sharedData = {};

function renderCIFSummary() {
  const numCIF = parseInt(document.querySelector("[name='num.cifs']").value);
  const diffCIFProp = document.querySelector("[name='diff.cif.fraction']").value;
  const numDiffCIF = Math.ceil(numCIF * diffCIFProp);
  const numNonDiffCIF = numCIF - numDiffCIF;

  const isDiscrete = document.querySelector("[name=_cellPop]:checked").value === "discrete";
  const ncells = parseInt(document.querySelector("[name='num.cells']").value);
  const cifMean = parseFloat(document.querySelector("[name='cif.center']").value);
  const cifSD = parseFloat(document.querySelector("[name='cif.sigma']").value);
  const treeName = document.querySelector("[name='tree']:checked").value;
  const givMean = parseFloat(document.querySelector("[name='giv.mean']").value);
  const givSD = parseFloat(document.querySelector("[name='giv.sd']").value);
  const givProb = parseFloat(document.querySelector("[name='giv.prob']").value);
  const nclus = parseInt(document.querySelector("[name='_numCluster']").value);

  const grnEnabled = document.querySelector("[name=_grnEnabled]:checked").value === "on";
  let grnNumGenes = 0;
  let grnNumReg = 0;
  if (grnEnabled && !sharedData.GRNInfo.error) {
    grnNumReg = sharedData.GRNInfo.nregulators;
    grnNumGenes = sharedData.GRNInfo.ngenes;
  }
  const ngenes = parseInt(document.querySelector(`[name='num.genes${grnEnabled ? "" : "2"}']`).value);
  const keys = ["nonDiff", "diff", "grn"];
  const nCIF = { nonDiff: numNonDiffCIF, diff: numDiffCIF, grn: grnNumReg };
  const totalCIF = numCIF + grnNumReg;
  const totalWidth = 160;
  const totalHeight = 180;
  const size = {};
  let cumX = 1;
  for (const k of keys) {
    const w = (nCIF[k] / totalCIF) * totalWidth;
    size[k] = { w, x: cumX };
    cumX += size[k].w;
    size[k].labelX = size[k].x + size[k].w / 2;
  }
  const grnGenesWidth = (grnNumGenes / ngenes) * totalHeight;
  size.grnGenes = {
    w: grnGenesWidth,
    labelX: grnGenesWidth / 2,
    remLabelX: grnGenesWidth + (totalHeight - grnGenesWidth) / 2,
  };

  const cifTitleGRN = grnEnabled ? ` + ${grnNumReg} (GRN)` : "";
  const givTitleGRN = grnEnabled ? ` (${grnNumGenes} in GRN)` : "";

  const svg = `
  <svg xmlns="http://www.w3.org/2000/svg" id="cif-svg" width="960" height="220">
    <g transform="translate(16,1)">
      <text x="-4" y="110" font-size="12" text-anchor="middle" transform="rotate(-90 -4,110)">Cells: ${ncells}</text>
      <g transform="translate(0,16)">
          <text x="80" y="-4" font-size="12" text-anchor="middle">
          <tspan font-weight="bold">CIF</tspan>: ${numCIF}${cifTitleGRN}
          </text>
          <rect x="0" y="0" width="162" height="182" fill="white" stroke="black" />
          <rect x="${size.nonDiff.x}" y="1" width="${size.nonDiff.w}" height="180" fill="rgb(207, 226, 255)"
              stroke="black" />
          <text x="${size.nonDiff.labelX}" y="194" text-anchor="middle" font-size="12">${nCIF.nonDiff}</text>
          <rect x="${size.diff.x}" y="1" width="${size.diff.w}" height="180" fill="rgb(209, 231, 221)"
              stroke="black" />
          <text x="${size.diff.labelX}" y="194" text-anchor="middle" font-size="12">${nCIF.diff}</text>
          ${
            grnEnabled
              ? `
          <rect x="${size.grn.x}" y="1" width="${size.grn.w}" height="180" fill="rgb(248, 215, 218)"
              stroke="black" />
          <text x="${size.grn.labelX}" y="194" text-anchor="middle" font-size="12">${nCIF.grn}</text>
          `
              : ""
          }
      </g>
      <g transform="translate(220,26)">
        <text x="-32" y="90" font-size="24" text-anchor="middle">╳</text>
        <text x="90" y="-4" font-size="12" text-anchor="middle">Genes: ${ngenes}${givTitleGRN}</text>
        <text x="-4" y="80" font-size="12" text-anchor="middle" transform="rotate(-90 -4,80)">
        <tspan font-weight="bold">GIV</tspan>: ${totalCIF}
        </text>
        <rect x="0" y="0" width="182" height="162" fill="white" stroke="black" />
        <rect y="${size.nonDiff.x}" x="1" height="${size.nonDiff.w}" width="180" fill="rgb(226, 217, 255)"
            stroke="black" />
        <rect y="${size.diff.x}" x="1" height="${size.diff.w}" width="180" fill="rgb(226, 217, 255)"
            stroke="black" />
        ${
          grnEnabled
            ? `
        <rect y="${size.grn.x}" x="1" height="${size.grn.w}" width="180" fill="rgb(226, 217, 255)"
            stroke="black" />
        <rect y="${size.grn.x + 1}" x="1" height="${size.grn.w - 2}"
          width="${size.grnGenes.w}" fill="#e2e2e2"/>
        <text x="${size.grnGenes.labelX}" y="${size.grn.labelX + 6}"
              text-anchor="middle" font-size="12">GRN</text>
        <rect y="1" x="1" height="160" width="${size.grnGenes.w}" fill="rgba(0,0,0,.05)" stroke="black" stroke-dasharray="2,2" />
        <text x="${size.grnGenes.labelX}" y="174" text-anchor="middle" font-size="12">${grnNumGenes}</text>
        <text x="${size.grnGenes.remLabelX}" y="174" text-anchor="middle" font-size="12">${ngenes - grnNumGenes}</text>
        `
            : ""
        }
      </g>
      <g transform="translate(460,16)">
        <text x="-32" y="100" font-size="24" text-anchor="middle">=</text>
        <text x="90" y="-4" font-size="12" text-anchor="middle">Genes: ${ngenes}</text>
        <text x="-4" y="90" font-size="12" text-anchor="middle" transform="rotate(-90 -4,90)">Cells: ${ncells}</text>
        <rect x="0" y="0" width="182" height="182" fill="white" stroke="black" stroke-width="2" />
        <text x="90" y="96" font-size="12" text-anchor="middle">Kinetic Parameters</text>
      </g>
      <g transform="translate(700,16)">
        <line x1="-24" y1="0" x2="-24" y2="180" stroke="black" stroke-dasharray="4,4" />
        <g>
          <rect x="0" y="8" width="20" height="20" fill="rgb(207, 226, 255)" stroke="black" />
          <text x="0" y="0" font-size="12" font-weight="bold">
          non-diff-CIF (${numNonDiffCIF}) <tspan fill="rgb(180,180,180)">Cellular Heterogenity</tspan>
          </text>
          <text x="24" y="22" font-size="14">~ <tspan fill="#0d6efd">N(${cifMean}, ${cifSD})</tspan></text>
        </g>
        <g transform="translate(0,50)">
          <rect x="0" y="8" width="20" height="20" fill="rgb(209, 231, 221)" stroke="black" />
          <text x="0" y="0" font-size="12" font-weight="bold">
          diff-CIF (${numDiffCIF}) <tspan fill="rgb(180,180,180)">Cell Population</tspan>
          </text>
          <text x="24" y="22" font-size="14">
          ${
            isDiscrete
              ? `~ sampled from <tspan fill="#fd7e14">${nclus} clusters</tspan>`
              : `~ sampled along <tspan fill="#fd7e14">tree: ${treeName}</tspan>`
          }
          </text>
        </g>
        <g transform="translate(0,100)">
          <rect x="0" y="8" width="20" height="20" fill="rgb(248, 215, 218)" stroke="black" />
          <text x="0" y="0" font-size="12" font-weight="bold">
          grn-CIF (${grnNumReg}) <tspan fill="rgb(180,180,180)">TF Expression</tspan>
          </text>
          <text x="24" y="22" font-size="14">~ populated during simulation</text>
        </g>
        <g transform="translate(0,150)">
          <rect x="0" y="8" width="20" height="20" fill="rgb(226, 217, 255)" stroke="black" />
          <text x="0" y="0" font-size="12" font-weight="bold">
          GIV <tspan fill="rgb(180,180,180)">How gene is affected by CIF entries</tspan>
          </text>
          <text x="24" y="22" font-size="14"> ~ 0 w.p; <tspan fill="#0d6efd">${givProb},</tspan> <tspan fill="#0d6efd">N(${givMean}, ${givSD})</tspan> otherwise</text>
        </g>
      </g
    </g>
  </svg>
  `;

  document.getElementById("cif-svg-container").innerHTML = svg;
}

function renderATACSummary() {
  const numCIF = parseInt(document.querySelector("[name='num.cifs']").value);
  const ncells = parseInt(document.querySelector("[name='num.cells']").value);
  const ngenes = parseInt(document.querySelector("[name='num.genes']").value);
  const rivMean = parseFloat(document.querySelector("[name='giv.mean']").value);
  const rivSD = parseFloat(document.querySelector("[name='giv.sd']").value);
  const rivProb = parseFloat(document.querySelector("[name='giv.prob']").value);
  const atacEffect = parseFloat(document.querySelector("[name='atac.effect']").value);
  const regionDists = [0, 1, 2].map((i) => parseFloat(document.querySelector(`[name='_regionDist${i}']`).value));
  const cumsum = regionDists.map((d, i) => regionDists.slice(0, i + 1).reduce((a, b) => a + b, 0));
  const geneCumsum = cumsum.map((c) => Math.round(c * ngenes));
  const geneNums = geneCumsum.map((c, i) => (i === 0 ? c : c - geneCumsum[i - 1]));
  const widths = regionDists.map((d) => d * 100);
  const geneLabelX = widths.map((w, i) => (i === 0 ? w / 2 : cumsum[i - 1] * 100 + w / 2));
  console.log(widths, geneLabelX, geneNums);

  const svg = `
  <svg xmlns="http://www.w3.org/2000/svg" id="atac-svg" width="960" height="280">
    <defs>
      <marker
        id="triangle"
        viewBox="0 0 4 4"
        refX="1"
        refY="2"
        markerUnits="strokeWidth"
        markerWidth="4"
        markerHeight="4"
        orient="auto">
        <path d="M 0 0 L 4 2 L 0 4 z" fill="black" />
      </marker>
    </defs>
    <g transform="translate(16,1)">
      <g transform="translate(0,1)">
        <text x="-4" y="66" font-size="12" text-anchor="middle" transform="rotate(-90 -4,66)">Cells: ${ncells}</text>
        <g transform="translate(0,16)">
          <text x="40" y="-4" font-size="12" text-anchor="middle">
          <tspan font-weight="bold">CIF</tspan>: ${numCIF}
          </text>
          <rect x="0" y="0" width="82" height="102" fill="white" stroke="black" />
          <rect x="1" y="1" width="80" height="100" fill="#e2e2e2" stroke="black" />
          <text x="40" y="56" font-size="12" text-anchor="middle">CIF</text>
        </g>
      </g>
      <g transform="translate(140,26)">
        <text x="-32" y="50" font-size="24" text-anchor="middle">╳</text>
        <text x="50" y="-4" font-size="12" text-anchor="middle">Regions: ${ngenes * 3}</text>
        <text x="-4" y="46" font-size="12" text-anchor="middle" transform="rotate(-90 -4,46)">
        <tspan font-weight="bold">RIV</tspan>: ${numCIF}
        </text>
        <rect x="0" y="0" width="102" height="82" fill="white" stroke="black" />
        <rect x="1" y="1" width="100" height="80" fill="#ffe5d0" stroke="black" />
        <text x="50" y="46" font-size="12" text-anchor="middle">RIV</text>
        <text x="0" y="100" font-size="12" text-anchor="start"> = 0 w.p. <tspan fill="#0d6efd">${rivProb}</tspan></text>
        <text x="0" y="116" font-size="12" text-anchor="start"> ~ <tspan fill="#0d6efd">N(${rivMean}, ${rivSD})</tspan> otherwise</text>
      </g>
      <g transform="translate(340,16)">
        <line x1="-82" y1="52" x2="-30" y2="52" stroke="black" stroke-width="2" marker-end="url(#triangle)" />
        <text x="50" y="-4" font-size="12" text-anchor="middle">Regions: ${ngenes * 3}</text>
        <text x="-4" y="56" font-size="12" text-anchor="middle" transform="rotate(-90 -4,56)">Cells: ${ncells}</text>
        <rect x="0" y="0" width="102" height="102" fill="#cff4fc" stroke="black" stroke-width="2" />
        <text x="50" y="56" font-size="12" text-anchor="middle">scATAC-seq</text>
      </g>
      <g transform="translate(500,16)">
        <text x="-32" y="60" font-size="24" text-anchor="middle">╳</text>
        <text x="50" y="-4" font-size="12" text-anchor="middle">Genes: ${ngenes}</text>
        <text x="-4" y="56" font-size="12" text-anchor="middle" transform="rotate(-90 -4,56)">
        Regions: ${ngenes * 3}
        </text>
        <rect x="1" y="1" width="100" height="100" fill="#d2f4ea" stroke="black" stroke-width="2" />
        <text x="50" y="56" font-size="12" text-anchor="middle">Region-Gene</text>
        <text x="${geneLabelX[0]}" y="114" font-size="12" text-anchor="middle">${geneNums[0]}</text>
        <text x="${geneLabelX[1]}" y="114" font-size="12" text-anchor="middle">${geneNums[1]}</text>
        <text x="${geneLabelX[2]}" y="114" font-size="12" text-anchor="middle">${geneNums[2]}</text>
        <line x1="${widths[0]}" y1="100" x2="${widths[0]}" y2="106" stroke="black" />
        <line x1="${widths[1] + widths[0]}" y1="100" x2="${widths[1] + widths[0]}" y2="106" stroke="black" />
        <g transform="translate(0,146)">
          <text x="-4" y="0" font-size="12" text-anchor="end">Controlled by 0 regions</text>
          <text x="-4" y="16" font-size="12" text-anchor="end">Controlled by 1 regions</text>
          <text x="-4" y="32" font-size="12" text-anchor="end">Controlled by 2 regions</text>
          <g transform="translate(0,-4)">
            <path d="M 0 0 L ${geneLabelX[0]} 0 L ${geneLabelX[0]} -22" stroke="black" fill="none" />
            <path d="M 0 16 L ${geneLabelX[1]} 16 L ${geneLabelX[1]} -22" stroke="black" fill="none" />
            <path d="M 0 32 L ${geneLabelX[2]} 32 L ${geneLabelX[2]} -22" stroke="black" fill="none" />
          </g>
        </g>
      </g>
      <g transform="translate(660,16)">
        <text x="-32" y="60" font-size="24" text-anchor="middle">=</text>
        <text x="50" y="-4" font-size="12" text-anchor="middle">Genes: ${ngenes}</text>
        <text x="-4" y="56" font-size="12" text-anchor="middle" transform="rotate(-90 -4,56)">
        Cells: ${ncells}
        </text>
        <rect x="1" y="1" width="100" height="100" fill="#f7d6e6" stroke="black" stroke-width="2" />
        <g transform="translate(0,56)">
          <text x="50" y="-32" font-size="12" text-anchor="middle">Kinetic</text>
          <text x="50" y="-16" font-size="12" text-anchor="middle">Parameter</text>
          <text x="50" y="0" font-size="12" font-weight="bold" text-anchor="middle">k_on</text>
          <text x="50" y="16" font-size="12" text-anchor="middle">(from ATAC)</text>
          <text x="50" y="32" font-size="12" text-anchor="middle">with effect <tspan fill="#0d6efd">${atacEffect}</tspan></text>
        </g>
      </g>
      <g transform="translate(660,166)">
        <text x="50" y="-20" font-size="24" text-anchor="middle">+</text>
        <text x="50" y="-4" font-size="12" text-anchor="middle">Genes: ${ngenes}</text>
        <text x="-4" y="56" font-size="12" text-anchor="middle" transform="rotate(-90 -4,56)">
        Cells: ${ncells}
        </text>
        <rect x="1" y="1" width="100" height="100" fill="#e2e2e2" stroke="black" stroke-width="2" />
        <line x1="-600" y1="52" x2="-20" y2="52" stroke="black" stroke-width="2" marker-end="url(#triangle)" stroke-dasharray="4,4" />
        <text x="-606" y="56" font-size="12" text-anchor="end">CIF/GIV</text>
        <g transform="translate(0,56)">
          <text x="50" y="-24" font-size="12" text-anchor="middle">Kinetic</text>
          <text x="50" y="-8" font-size="12" text-anchor="middle">Parameter</text>
          <text x="50" y="8" font-size="12" font-weight="bold" text-anchor="middle">k_on</text>
          <text x="50" y="24" font-size="12" text-anchor="middle">(from CIF/GIV)</text>
        </g>
      </g>
      <g transform="translate(840,90)">
        <line x1="-62" y1="52" x2="-30" y2="52" stroke="black" stroke-width="2" marker-end="url(#triangle)" />
        <text x="50" y="-4" font-size="12" text-anchor="middle">Genes: ${ngenes}</text>
        <text x="-4" y="56" font-size="12" text-anchor="middle" transform="rotate(-90 -4,56)">
        Cells: ${ncells}
        </text>
        <rect x="1" y="1" width="100" height="100" fill="#e2e2e2" stroke="black" stroke-width="2" />
        <g transform="translate(0,56)">
          <text x="50" y="-24" font-size="12" text-anchor="middle">Kinetic</text>
          <text x="50" y="-8" font-size="12" text-anchor="middle">Parameter</text>
          <text x="50" y="8" font-size="12" font-weight="bold" text-anchor="middle">k_on</text>
          <text x="50" y="24" font-size="12" text-anchor="middle">Actually used</text>
        </g>
      </g>
    </g>
  </svg>`;
  document.getElementById("atac-svg-container").innerHTML = svg;
}

function updateTreeSummary() {
  if (sharedData.TreeInfo.error) {
    document.getElementById("tree_summary").innerHTML = sharedData.TreeInfo.error;
    return;
  }
  const ntips = Array.isArray(sharedData.TreeInfo.tips) ? sharedData.TreeInfo.tips.length : 1;
  const nnodes = sharedData.TreeInfo.internal.length;
  const nedges = sharedData.TreeInfo.edges.length;
  let text = `Tree with ${nnodes} nodes, ${nedges} edges, and ${ntips} leaves.`;
  if (ntips < 2) {
    text += "<br>(tree with less than 2 tips may not be plotted)";
  }
  document.getElementById("tree_summary").innerHTML = text;
}

function updateSpatialPreview() {
  const grid = sharedData.Grid;
  const svgSize = 600;
  const dotSize = svgSize / grid.size;
  const ctypes = Array.from(new Set(grid.final_types));
  const colors = ctypes.map((_, i) => colorList[i % colorList.length]);
  const colorMap = ctypes.reduce((acc, c, i) => ({ ...acc, [c]: colors[i] }), {});
  const icells = grid.final_types.map((_, i) => i);
  const svg = `
  <svg xmlns="http://www.w3.org/2000/svg" id="spatial-svg" width="${svgSize}" height="${svgSize}">
    <g transform="translate(-${dotSize / 2},-${dotSize / 2})">
    ${icells.map((i) => {
      const x = grid.locs[i][0] * dotSize;
      const y = grid.locs[i][1] * dotSize;
      return `<circle cx="${x}" cy="${y}" r="${dotSize / 2}" fill="${colorMap[grid.final_types[i]]}" />`;
    })}
    </g>
  </svg>
  `;
  document.getElementById("spatial-svg-container").innerHTML = svg;
}

function setupSpatialPreviewButton() {
  let idx = 0;
  const btn = document.getElementById("submit_spatial");
  const msgDiv = document.getElementById("submit_spatial_msg");

  btn.addEventListener("click", function () {
    const gridSize = parseInt(document.querySelector("[name='sp.grid.size']").value);
    const ncell = parseInt(document.querySelector("[name='num.cells']").value);
    let layout = document.querySelector("[name='_spLayout']").value;
    const sameTypeProb = parseFloat(document.querySelector("[name='sp.same.type.prob']").value);
    const stepSize = parseFloat(document.querySelector("[name='sp.step.size']").value);
    if (layout === "islands") {
      const spIslands = document.querySelector("[name='_spIslands']").value;
      layout = `${layout}:${spIslands}`;
    } else if (layout === "normal") {
      layout = "enhanced";
    }
    if (ncell > gridSize ** 2) {
      msgDiv.innerText = "Number of cells exceeds grid size";
      msgDiv.classList.add("text-danger");
      return;
    }
    msgDiv.classList.remove("text-danger");
    Shiny.setInputValue("submit_spatial", {
      idx,
      ncell,
      gridSize,
      layout,
      sameTypeProb,
      stepSize,
    });
    btn.disabled = true;
    msgDiv.innerText = "Generating...";
    idx += 1;
  });

  Shiny.addCustomMessageHandler("Grid", (grid) => {
    sharedData.Grid = grid;
    btn.disabled = false;
    msgDiv.innerText = "";
    console.log(grid);
    updateSpatialPreview();
  });
}

function init() {
  Shiny.addCustomMessageHandler("RObjects", (objList) => {
    console.log(objList);
    const grnSelect = document.getElementById("R_grnSelect");
    const cciSelect = document.getElementById("R_cciSelect");
    [grnSelect, cciSelect].forEach((e, i) => {
      const name = i === 0 ? "GRN" : "_CCI";
      e.innerHTML = objList
        .flatMap((obj) => (obj.type === "data.frame" ? `<li><a class="dropdown-item" href="#">${obj.name}</a></li>` : []))
        .join("");
      e.querySelectorAll("a").forEach((a) =>
        a.addEventListener("click", function (e) {
          const newValue = e.target.innerText;
          document.querySelector(`[name=${name}]`).value = newValue;
          Shiny.setInputValue(name, newValue);
        })
      );
    });
  });

  Shiny.addCustomMessageHandler("GRNInfo", (info) => {
    sharedData.GRNInfo = info;
    renderCIFSummary();
    renderATACSummary();
  });

  Shiny.addCustomMessageHandler("TreeInfo", (info) => {
    sharedData.TreeInfo = info;
    updateTreeSummary();
  });

  function updateTreePreview() {
    let treeVal = document.querySelector("[name=tree]:checked").value;
    if (treeVal === "custom") {
      treeVal = document.querySelector("[name='_treeCustom']").value;
    }
    Shiny.setInputValue("tree", treeVal);
  }

  document.querySelectorAll("[name=tree]").forEach((el) => {
    el.addEventListener("change", updateTreePreview);
  });
  document.querySelector("[name='_treeCustom']").addEventListener("change", updateTreePreview);
  setTimeout(updateTreePreview, 1000);

  document.querySelectorAll("[name=_grnEnabled]").forEach((el) =>
    el.addEventListener("change", function (e) {
      const id = e.target.value === "on" ? "#grn-enable-tab" : "#grn-disable-tab";
      bootstrap.Tab.getOrCreateInstance(document.querySelector(id)).show();
    })
  );
  document.querySelectorAll("[name=_cellPop]").forEach((el) =>
    el.addEventListener("change", function (e) {
      const id = e.target.value === "continuous" ? "#pop-cont-tab" : "#pop-disc-tab";
      bootstrap.Tab.getOrCreateInstance(document.querySelector(id)).show();
    })
  );
  document.querySelectorAll("[name=velocity]").forEach((el) =>
    el.addEventListener("change", function (e) {
      const id = e.target.value === "on" ? "#velo-enable-tab" : "#velo-disable-tab";
      bootstrap.Tab.getOrCreateInstance(document.querySelector(id)).show();
    })
  );
  document.querySelector("[name=_spLayout]").addEventListener("change", function (e) {
    const id = `#sp-${e.target.value}-tab`;
    console.log(id);
    bootstrap.Tab.getOrCreateInstance(document.querySelector(id)).show();
  });

  observe(
    ["num.cifs", "diff.cif.fraction", "_cellPop", "num.cells", "cif.center", "cif.sigma", "tree", "giv.mean", "giv.sd", "giv.prob"],
    renderCIFSummary
  );
  observe(
    ["num.cifs", "num.genes", "num.cells", "riv.mean", "riv.sd", "riv.prob", "atac.effect", "_regionDist0", "_regionDist1", "_regionDist2"],
    renderATACSummary
  );

  setupSpatialPreviewButton();
}

function observe(names, func) {
  const els = names.flatMap((name) => Array.from(document.querySelectorAll(`[name='${name}']`)));
  els.forEach((el) => el.addEventListener("change", func));
}

document.addEventListener("DOMContentLoaded", init);
