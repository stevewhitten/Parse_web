const window_size = 25;
const amino_acids = {
    A: { ppii_h: 0.37, beta: 0.77 },
    C: { ppii_h: 0.25, beta: 0.81 },
    D: { ppii_h: 0.3, beta: 1.41 },
    E: { ppii_h: 0.42, beta: 0.99 },
    F: { ppii_h: 0.17, beta: 0.59 },
    G: { ppii_h: 0.13, beta: 1.64 },
    H: { ppii_h: 0.2, beta: 0.68 },
    I: { ppii_h: 0.39, beta: 0.51 },
    K: { ppii_h: 0.56, beta: 0.96 },
    L: { ppii_h: 0.24, beta: 0.58 },
    M: { ppii_h: 0.36, beta: 0.41 },
    N: { ppii_h: 0.27, beta: 1.28 },
    P: { ppii_h: 1.0, beta: 1.91 },
    Q: { ppii_h: 0.53, beta: 0.98 },
    R: { ppii_h: 0.38, beta: 0.88 },
    S: { ppii_h: 0.24, beta: 1.32 },
    T: { ppii_h: 0.32, beta: 1.04 },
    V: { ppii_h: 0.39, beta: 0.47 },
    W: { ppii_h: 0.25, beta: 0.76 },
    Y: { ppii_h: 0.25, beta: 1.05 },
};
const amino_acid_keys = Object.keys(amino_acids);

function getSummary(input) {
    var count = {
        A: 0,
        C: 0,
        D: 0,
        E: 0,
        F: 0,
        G: 0,
        H: 0,
        I: 0,
        K: 0,
        L: 0,
        M: 0,
        N: 0,
        P: 0,
        Q: 0,
        R: 0,
        S: 0,
        T: 0,
        V: 0,
        W: 0,
        Y: 0,
    };

    for (let i = 0; i < input.length; i++) {
        count[input[i]]++;
    }

    var net_charge = Math.abs(count.D + count.E - count.K - count.R);

    var sum_ppii = 0.0;

    amino_acid_keys.forEach((x) => {
        sum_ppii += count[x] * amino_acids[x].ppii_h;
    });

    var fppii = sum_ppii / input.length;
    var v_exponent = 0.503 - 0.11 * Math.log(1.0 - fppii);

    var rh = 2.16 * input.length ** v_exponent + 0.26 * net_charge - 0.29 * input.length ** 0.5;
    var v_flory = Math.log(rh / 2.16) / Math.log(input.length);

    var turn = 0.0;
    amino_acid_keys.forEach((x) => {
        turn += count[x] * amino_acids[x].beta;
    });

    turn = turn / input.length;

    var result = {
        nu: v_flory,
        turn: turn,
        r_model: turn / v_flory,
    };
    return result;
}

function getResidues(input) {
    var residues = [];
    for (let i = 0; i < input.length; i++) {
        residues.push({
            amino_acid: input[i],
            domain: null,
            r_model: null,
        });
    }

    for (let i = 0; i < input.length; i++) {
        if (i < input.length - window_size + 1) {
            var middle_position = i + Math.trunc(window_size / 2);
            var count = {
                A: 0,
                C: 0,
                D: 0,
                E: 0,
                F: 0,
                G: 0,
                H: 0,
                I: 0,
                K: 0,
                L: 0,
                M: 0,
                N: 0,
                P: 0,
                Q: 0,
                R: 0,
                S: 0,
                T: 0,
                V: 0,
                W: 0,
                Y: 0,
            };

            for (let j = i; j < i + window_size; j++) {
                count[input[j]]++;
            }

            var net_charge = Math.abs(count.D + count.E - count.K - count.R);

            var sum_ppii = 0.0;

            amino_acid_keys.forEach((x) => {
                sum_ppii += count[x] * amino_acids[x].ppii_h;
            });

            var fppii = sum_ppii / window_size;
            var v_exponent = 0.503 - 0.11 * Math.log(1.0 - fppii);

            var rh = 2.16 * (4 * window_size) ** v_exponent + 0.26 * 4 * net_charge - 0.29 * (4 * window_size) ** 0.5;
            var v_flory = Math.log(rh / 2.16) / Math.log(4 * window_size);

            var turn = 0.0;
            amino_acid_keys.forEach((x) => {
                turn += count[x] * amino_acids[x].beta;
            });

            turn = turn / window_size;

            var slope_1 = 0.019 / 0.082;
            var intercept_1 = 0.539 - slope_1 * 1.062;
            var v_line_1 = slope_1 * turn + intercept_1;

            var slope_2 = -0.019 / 0.082;
            var intercept_2 = 0.539 - slope_2 * 1.062;
            var v_line_2 = slope_2 * turn + intercept_2;

            residues[middle_position].r_model = turn / v_flory;
            if (v_flory >= 0.558) {
                residues[middle_position].domain = 'D';
            } else if (turn > 1.144 && v_flory < 0.558) {
                residues[middle_position].domain = 'P';
            } else if (turn > 1.062 && v_flory < 0.539) {
                residues[middle_position].domain = 'P';
            } else if (turn > 1.062 && v_flory < v_line_1) {
                residues[middle_position].domain = 'P';
            } else if (turn > 1.062 && v_flory >= v_line_1) {
                residues[middle_position].domain = 'D';
            } else if (turn <= 0.98 && v_flory < 0.558) {
                residues[middle_position].domain = 'F';
            } else if (turn <= 1.062 && v_flory < 0.539) {
                residues[middle_position].domain = 'F';
            } else if (turn <= 1.062 && v_flory < v_line_2) {
                residues[middle_position].domain = 'F';
            } else if (turn <= 1.062 && v_flory >= v_line_2) {
                residues[middle_position].domain = 'D';
            }
        }
    }

    for (let i = 0; i < Math.trunc(window_size / 2); i++) {
        residues[i].domain = residues[Math.trunc(window_size / 2)].domain;
        residues[i].r_model = residues[Math.trunc(window_size / 2)].r_model;

        residues[input.length - i - 1].domain = residues[input.length - Math.trunc(window_size / 2) - 1].domain;
        residues[input.length - i - 1].r_model = residues[input.length - Math.trunc(window_size / 2) - 1].r_model;
    }

    return residues;
}

function getDomainRegion(residues, regionLength, cutoff, domain) {
    var i = 0;
    var regions = [];

    while (i + regionLength <= residues.length) {
        var j = i;
        var found_region = false;
        var count_p = 0;
        var count_w = 0;

        while (j < residues.length) {
            count_w++;
            if (residues[j].domain == domain) {
                count_p++;
            }

            if (count_w >= regionLength) {
                var percent_p = count_p / count_w;
                if (percent_p >= cutoff) {
                    found_region = true;
                } else {
                    break;
                }
            }

            j++;
        }

        if (found_region) {
            regions.push({
                domain: domain,
                start: i,
                end: j - 1,
            });
            i = j;
        } else {
            i++;
        }
    }

    return regions;
}

function getRegions(residues, regionLength, cutoff) {
    var regions = [];
    var regionList;

    regionList = getDomainRegion(residues, regionLength, cutoff, 'P');
    for (const x of regionList) {
        regions.push(x);
    }
    regionList = getDomainRegion(residues, regionLength, cutoff, 'D');
    for (const x of regionList) {
        regions.push(x);
    }
    regionList = getDomainRegion(residues, regionLength, cutoff, 'F');
    for (const x of regionList) {
        regions.push(x);
    }

    regions.sort((a, b) => {
        return a.start - b.start;
    });

    if (regions.length > 1) {
        for (let i = 0; i < regions.length - 1; i++) {
            if (regions[i].end >= regions[i + 1].start) {
                regions[i].end = regions[i + 1].start + Math.round(((1.0 - cutoff) * regionLength) / 2.0);
                regions[i + 1].start = regions[i].end + 1;
            }
        }
    }

    return regions;
}

function parse(sequence, regionLength, cutoff) {
    var result = {
        error: false,
        message: null,
        data: null,
    };

    if (!sequence) {
        result.error = true;
        result.message = `No sequence specified`;
        return result;
    }

    sequence = sequence.replace(/\s+/g, '');
    sequence = sequence.toUpperCase();

    sequence = Array.from(sequence);

    if (sequence.length < 25) {
        result.error = true;
        result.message = `Sequence length should be at least 25 characters long`;
        return result;
    } else if (sequence.length > 10000) {
        result.error = true;
        result.message = `Sequence length should be at most 10000 characters long`;
        return result;
    }

    var summary = getSummary(sequence);
    var residues = getResidues(sequence);
    var regions = getRegions(residues, regionLength, cutoff);

    var data = {
        sequence: sequence,
        summary: summary,
        residues: residues,
        regions: regions,
    };

    result.data = data;

    return result;
}
