const inputSequenceEl = document.getElementById('input-sequence');
const inputFileEl = document.getElementById('input-file');
const inputRegionLengthEl = document.getElementById('input-region-length');
const inputCutoffEl = document.getElementById('input-cutoff');
const submitEl = document.getElementById('submit');
const clearEl = document.getElementById('clear');
const resetEl = document.getElementById('reset');
const radioButtons = document.querySelectorAll("input[name='display']");

const downloadRegionsEl = document.getElementById('download-regions-csv');
const downloadResiduesEl = document.getElementById('download-residues-csv');

const proteinsEl = document.getElementById('proteins');
const proteinInformationEl = document.getElementById('protein-information');
const setInformationEl = document.getElementById('set-information');

const outputSequenceLengthEl = document.getElementById('output-sequence-length');
const outputSequenceNuModelEl = document.getElementById('output-sequence-nu-model');
const outputSequenceBetaTurnEl = document.getElementById('output-sequence-beta-turn');
const outputSequenceRModelEl = document.getElementById('output-sequence-r-model');

const proteinChartEl = document.getElementById('protein-chart');
const setChartEl = document.getElementById('set-chart');

var selectedDisplay = 'protein';
const proteins = [];

const regionsTable = new Tabulator('#regions-table', {
    autoResize: true,
    maxHeight: '300px',
    layout: 'fitColumns',
    columns: [
        {
            title: 'Region number',
            field: 'index',
            sorter: 'number',
        },
        {
            title: 'Region classifier',
            field: 'domain',
            sorter: 'string',
        },
        {
            title: 'First residue of region',
            field: 'start',
            sorter: 'number',
        },
        {
            title: 'Last residue of region',
            field: 'end',
            sorter: 'number',
        },
        {
            title: 'Region length',
            field: 'length',
            sorter: 'number',
        },
    ],
});
const residuesTable = new Tabulator('#residues-table', {
    autoResize: true,
    maxHeight: '300px',
    layout: 'fitColumns',
    columns: [
        {
            title: 'Residue number',
            field: 'index',
            sorter: 'number',
        },
        {
            title: 'Amino Acid type',
            field: 'amino_acid',
            sorter: 'string',
        },
        {
            title: 'Residue label',
            field: 'domain',
            sorter: 'number',
        },
        {
            title: '<i>r<sub>model</sub></i>',
            field: 'r_model',
            sorter: 'number',
        },
    ],
});

const proteinChartLayout = {
    showlegend: true,
    xaxis: {
        title: {
            text: 'Residue Number',
        },
    },
    yaxis: {
        title: {
            text: 'Residue-level  <i>r<sub>model</sub></i>',
        },
    },
    shapes: [
        {
            type: 'line',
            xref: 'paper',
            yref: 'y',
            x0: 0,
            x1: 1,
            y0: 1.58,
            y1: 1.58,
            line: {
                dash: 'dash',
            },
        },
        {
            type: 'line',
            xref: 'paper',
            yref: 'y',
            x0: 0,
            x1: 1,
            y0: 1.88,
            y1: 1.88,
        },
        {
            type: 'line',
            xref: 'paper',
            yref: 'y',
            x0: 0,
            x1: 1,
            y0: 1.28,
            y1: 1.28,
        },
    ],
};
const proteinChartData = [
    {
        showlegend: false,
        mode: 'text',
        textposition: 'top right',
        text: ['huProteome mean = 1.58 ± 0.3'],
        hoverinfo: 'skip',
        x: [0],
        y: [1.58],
    },
];

const setChartLayout = {
    showlegend: false,
    xaxis: {
        title: {
            text: 'Predicted PS region length',
        },
    },
    yaxis: {
        title: {
            text: '% of set',
        },
        range: [0, 100],
    },
};

function resetSummary() {
    outputSequenceLengthEl.innerText = '';
    outputSequenceNuModelEl.innerText = '';
    outputSequenceBetaTurnEl.innerText = '';
    outputSequenceRModelEl.innerText = '';
}
function resetPlots() {
    Plotly.newPlot(proteinChartEl, [], proteinChartLayout);
    Plotly.newPlot(setChartEl, [], proteinChartLayout);
}
function resetTables() {
    regionsTable.setData([]);
    residuesTable.setData([]);
}

function updateDisplay(event) {
    if (event.target.value == 'set') {
        proteinInformationEl.setAttribute('hidden', true);
        setInformationEl.removeAttribute('hidden');
        Plotly.newPlot(setChartEl, setChartEl.data, setChartEl.layout);
    } else {
        setInformationEl.setAttribute('hidden', true);
        proteinInformationEl.removeAttribute('hidden');
        Plotly.newPlot(proteinChartEl, proteinChartEl.data, proteinChartEl.layout);
    }
}

function updateProteinSelect() {
    proteinsEl.innerHTML = '';
    proteins.forEach((x, index) => {
        var option = document.createElement('option');
        option.value = index;
        option.text = x.id;
        proteinsEl.appendChild(option);
    });
    if (proteins.length == 1) {
        proteinsEl.setAttribute('hidden', true);
        proteinsEl.removeAttribute('hidden');
    }
}

function selectProtein(event) {
    const selectedProtein = proteins[event.target.value];

    outputSequenceLengthEl.innerText = selectedProtein.sequence.length;
    outputSequenceNuModelEl.innerText = selectedProtein.summary.nu.toFixed(3);
    outputSequenceBetaTurnEl.innerText = selectedProtein.summary.turn.toFixed(3);
    outputSequenceRModelEl.innerText = selectedProtein.summary.r_model.toFixed(3);

    var regionData = [];
    selectedProtein.regions.forEach((x, index) => {
        regionData.push({
            index: index + 1,
            domain: x.domain,
            start: x.start + 1,
            end: x.end + 1,
            length: (x.end - x.start) + 1,
        });
    });
    regionsTable.setData(regionData);

    var residueData = [];
    selectedProtein.residues.forEach((x, index) => {
        residueData.push({
            index: index + 1,
            amino_acid: x.amino_acid,
            domain: x.domain,
            r_model: x.r_model.toFixed(3),
        });
    });
    residuesTable.setData(residueData);

    const traces = {
        all: {
            name: 'All',
            type: 'scatter',
            mode: 'lines',
            hoverinfo: 'none',
            showlegend: false,
            line: {
                color: '#c0c0c0',
            },
            x: [],
            y: [],
        },
        P: {
            name: 'P labeled residue',
            type: 'scatter',
            mode: 'lines',
            hoverinfo: 'x+y',
            line: {
                color: '#5075B0',
            },
            x: [],
            y: [],
        },
        D: {
            name: 'D labeled residue',
            type: 'scatter',
            mode: 'lines',
            hoverinfo: 'x+y',
            line: {
                color: '#A92217',
            },
            x: [],
            y: [],
        },
        F: {
            name: 'F labeled residue',
            type: 'scatter',
            mode: 'lines',
            hoverinfo: 'x+y',
            line: {
                color: '#000000',
            },
            x: [],
            y: [],
        },
    };

    selectedProtein.residues.forEach((x, index) => {
        traces.all.x.push(index);
        traces.all.y.push(x.r_model);
        traces.P.x.push(index);
        traces.P.y.push(x.domain == 'P' ? x.r_model : null);
        traces.D.x.push(index);
        traces.D.y.push(x.domain == 'D' ? x.r_model : null);
        traces.F.x.push(index);
        traces.F.y.push(x.domain == 'F' ? x.r_model : null);
    });

    const data = [...proteinChartData, traces.all, traces.P, traces.D, traces.F];

    Plotly.newPlot(proteinChartEl, data, proteinChartLayout);
}

function processProteinSequence(sequence, id) {
    var regionLength = parseInt(inputRegionLengthEl.value);
    var cutoff = parseFloat(inputCutoffEl.value);
    cutoff /= 100;

    var parseResult = parse(sequence, regionLength, cutoff);

    if (parseResult.error) {
        console.log(parseResult.message);
        console.log(id);
    } else {
        proteins.push({
            id: id,
            sequence: parseResult.data.sequence,
            summary: parseResult.data.summary,
            regions: parseResult.data.regions,
            residues: parseResult.data.residues,
        });
    }
}

function parseProteinSet(cutoff) {
    var cutoff = parseFloat(inputCutoffEl.value);
    cutoff /= 100;

    const setChartData = {
        //name: 'All',
        type: 'scatter',
        mode: 'lines',
        showlegend: false,
        line: {
            color: '#000000',
        },
        x: [],
        y: [],
    };

    if (proteins.length > 1) {
        var maxSequenceLength = 0;
        var proteinPHistogram = [];
        proteins.forEach((x) => {
            if (x.sequence.length > maxSequenceLength) {
                maxSequenceLength = x.sequence.length;
            }

            var longestRegion = 0;
            var i = 0;
            while (i < x.residues.length) {
                var j = i;
                var found_region = false;
                var count_p = 0;
                var count_w = 0;

                while (j < x.residues.length) {
                    count_w++;
                    if (x.residues[j].domain == 'P') {
                        count_p++;
                    }

                    var percent_p = count_p / count_w;
                    if (percent_p >= cutoff) {
                        found_region = true;
                    } else {
                        break;
                    }

                    j++;
                }

                if (found_region) {
                    var length = j - 1 - i;
                    if (length > longestRegion) {
                        longestRegion = length;
                    }
                    i = j;
                } else {
                    i++;
                }
            }

            if (!proteinPHistogram[longestRegion]) {
                proteinPHistogram[longestRegion] = 0;
            }
            proteinPHistogram[longestRegion]++;
        });

        for (let i = 0; i < maxSequenceLength; i++) {
            if (!proteinPHistogram[i]) {
                proteinPHistogram[i] = 0;
            }
        }

        for (let i = 0; i < maxSequenceLength; i++) {
            for (let j = i + 1; j < maxSequenceLength; j++) {
                proteinPHistogram[i] += proteinPHistogram[j];
            }
        }

        console.log(proteinPHistogram);
        for (let i = 0; i < 150; i++) {
            setChartData.x[i] = i;
            if (proteinPHistogram[i]) {
                setChartData.y[i] = (proteinPHistogram[i] / proteins.length) * 100;
            } else {
                setChartData.y[i] = null;
            }
        }
    }

    Plotly.newPlot(setChartEl, [setChartData], setChartLayout);
}

function parseFile(data) {
    var proteinsText = data.split(/^(?:>)/gm);
    proteinsText.forEach((x) => {
        if (x.length > 0) {
            proteinLines = x.split(/(?:\r\n|\n)+/g);
            proteinId = proteinLines[0];
            sequenceLines = proteinLines.splice(1);
            sequence = sequenceLines.join('');
            processProteinSequence(sequence, proteinId);
        }
    });
    return proteinsText.length;
}

function handleFile(file) {
    const fileReader = new FileReader();
    fileReader.addEventListener('load', (event) => {
        var data = fileReader.result;
        if (file.type == 'application/x-gzip') {
            try {
                data = pako.inflate(data);
                var textDecoder = new TextDecoder();
                data = textDecoder.decode(data);
            } catch (err) {
                console.log(err);
                alert('Could not decompress file');
                return;
            }
        }
        proteins.length = 0;
        var totalProteins = parseFile(data);
        if (proteins.length != totalProteins) {
            // Silently drop proteins that didnt parse
            // alert(`Could not parse ${totalProteins - proteins.length} sequences`);
        }
        console.log(`${proteins.length} proteins parsed`);
        parseProteinSet();
        updateProteinSelect();
        proteinsEl.dispatchEvent(new Event('change'));
    });
    if (file.type == 'application/x-gzip') {
        fileReader.readAsArrayBuffer(file);
    } else {
        fileReader.readAsText(file);
    }
}

function submit() {
    var sequence = inputSequenceEl.value;
    var files = inputFileEl.files;

    if (!sequence && files.length > 0) {
        const file = files[0];
        handleFile(file);
    } else {
        proteins.length = 0;
        processProteinSequence(sequence, '');
        if (proteins.length != 1) {
            alert('Could not parse sequence');
        } else {
            parseProteinSet();
            updateProteinSelect();
            proteinsEl.dispatchEvent(new Event('change'));
        }
    }
}
function clear() {
    console.log('clear');
    inputSequenceEl.value = '';
}
function reset() {
    clear();
    resetSummary();
    resetPlots();
    resetTables();
}

reset();

submitEl.addEventListener('click', submit);
clearEl.addEventListener('click', clear);
resetEl.addEventListener('click', reset);
proteinsEl.addEventListener('change', selectProtein);
radioButtons.forEach((x) => {
    x.addEventListener('change', updateDisplay);
});
downloadRegionsEl.addEventListener('click', _ => {
    regionsTable.download("csv", "regions.csv");
});
downloadResiduesEl.addEventListener('click', _ => {
    residuesTable.download("csv", "residues.csv");
});
