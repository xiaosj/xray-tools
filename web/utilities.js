// Generate an evenly spaced list, inclusive 
function linspace(start, end, step_num=100) {
    if(step_num > 2) {
        var list = [];
        var nstep = Math.floor(step_num);
        var step_value = (end - start) / (nstep - 1)
        for(var i = 0; i < nstep; i++) {
            list.push(start + step_value * i);
        }
        return list;
    } else {
        return [];
    }
}

// Generate a log-evenly spaced list, inclusive 
function linspace(start, end, step_num=100) {
    if(step_num > 2 && start > 0 && end > 0) {
        var log_start = Math.log(start);
        var log_end = Math.log(end);
        var nstep = Math.floor(step_num);
        var list = [];
        var step_value = (log_end - log_start) / (nstep - 1)
        for(var i = 0; i < nstep; i++) {
            list.push(Math.exp(log_start + step_value * i));
        }
        return list;
    } else {
        return [];
    }
}

// Interpolation functions for different linear / log choices
function lininterp(x1, x2, y1, y2, x0) {
    return y1 + (y2 - y1) / (x2 - x1) * (x0 - x1);
}

function loginterp(x1, x2, y1, y2, x0) {
    return Math.exp(lininterp(Math.log(x1), Math.log(x2), Math.log(y1), Math.log(y2), Math.log(x0)));
}

function linloginterp(x1, x2, y1, y2, x0) {
    return Math.exp(lininterp(x1, x2, Math.log(y1), Math.log(y2), x0));
}

function loglininterp(x1, x2, y1, y2, x0) {
    return lininterp(Math.log(x1), Math.log(x2), y1, y2, Math.log(x0));
}

// Interpolate from x / y list
function interp1d(x, y, x0, logx = false, logy = false) {
    var n = x.length;
    if(x0 < x[0] || x0 > x[n-1]) {
        return NaN;
    } else {
        // Find index
        var idx = Math.floor(n / 2);
        var idx0 = 0, idx1 = n - 1;
        var find = false;
        while (!find) {
            if (x[idx] <= x0 && x0 < x[idx+1]) {
                find = true;
            } else if (x0 < x[idx]) {
                idx1 = idx;
                idx = Math.floor((idx + idx0) / 2);
            } else {
                idx0 = idx;
                idx = Math.floor((idx + idx1) / 2);
            }
        }
        // Interpolate according to the linear / log choices
        var x1 = x[idx], x2 = x[idx + 1];
        var y1 = y[idx], y2 = y[idx + 1];
        if(!logx) {
            if(!logy) {
                lininterp(x1, y1, x2, y2, x0);
            } else {
                linloginterp(x1, y1, x2, y2, x0);
            }
        } else {
            if(!logy) {
                loglininterp(x1, y1, x2, y2, x0);
            } else {
                loginterp(x1, y1, x2, y2, x0);
            }
        }
    }
}
