// Generate an evenly spaced list, inclusive 
function linspace(start, end, step_num) {
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
function linspace(start, end, step_num) {
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
