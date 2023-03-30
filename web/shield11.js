// SHIELD11 input parameters
var Ebeam, NamTar, TarLen, TarRad, IattSW = false, FTDSW = 1;
var NamShl, AngShl, DisShl, ThkShl;
var Theta0, Theta1, Ntheta, minTheta = 5;
var AddCON, AddFE, AddPB, AddMIS;
var IDShld = 0, IDtar = 1;

var dose;
// 2D array to store all doses, the other dimension intialized in validate():
// 0: GNR;  1: MID; 2: HEN -- Neutrons
// 3: GamD; 4: GamI        -- Photons
// 5: total neutron; 6: total photon
// 7: total dose;    8: G/N ratio

// Parameters for web application
var powerUnit, units, oldunits, autoConv;
var eps;    // Electron per second
var doseUnit, dosenorm;   // Dose normalizaiton factor

var degree = Math.PI / 180.0, inch = 2.54, cm=1.0;
var lunitname, lunit = inch;

// Global variables for model drawing
var xminr, xmaxr, yminr, ymaxr; // x/y min/max in real world

function updateAll() {
  validate();
  calculate();
  updateDose();
  DrawDose();
}


//===========================================================
// UI functions
//===========================================================

// Convert between power and current
function powerModeChange() {
  var newPowerUnit = document.getElementById('powerMode').value;
  if(powerUnit != newPowerUnit) {
    if(newPowerUnit == "Power") {
      document.getElementById('powerunit').innerHTML = 'W';
      document.getElementById('power').value = (eps * Ebeam*1.e9 * 1.602e-19).toPrecision(3);
    }
    if(newPowerUnit == "Current") {
      document.getElementById('powerunit').innerHTML = 'mA';
      document.getElementById('power').value = (eps * 1.602e-19 * 1000).toPrecision(3);
    }
    powerUnit = newPowerUnit
  }
}

// Convert between SI and US units
function unitChange() {
  oldunits = units;
  if(document.getElementById('unit_SI').checked) {
    units = 'unit_SI';
    doseUnit = 'uSv/h';
    dosenorm = eps * 1e4 * 3600;
    lunit = cm;
    lunitname = 'cm';
  }
  else {
    units = 'unit_US';
    doseUnit = 'mrem/h';
    dosenorm = eps * 1e3 * 3600;
    lunit = inch;
    lunitname = 'inch';
  }

  // Update length unit
  var ltext = document.getElementsByClassName('lunit');
  for (var i = 0; i < ltext.length; i++) {  
    ltext[i].innerHTML = lunitname;
  };

  // Update dose unit
  var dtext = document.getElementsByClassName('dunit');
  for (var i = 0; i < dtext.length; i++) {
    dtext[i].innerHTML = doseUnit;
  };

  autoConv = document.getElementById('autoConv').checked;
  if(autoConv) {  // Covert length values according to the selected unit
    var linp = document.getElementsByClassName('inpL');
    var unitfactor;
    if(units=='unit_SI' && oldunits=='unit_US')
      unitfactor = inch/cm;
    else if(units=='unit_US' && oldunits=='unit_SI')
      unitfactor = cm/inch;
    else
      unitfactor = 1.0;
    for (var i = 0; i < linp.length; i++) {
      linp[i].value = (parseFloat(linp[i].value) * unitfactor).toFixed(2);
    };
  }
  updateAll();
}

// Calculate electrons per second
function GetEPS() {
  var inpv = parseFloat(document.getElementById('power').value);
  if(inpv <= 0) {
    inpv = 1
    document.getElementById('power').value = '1'
  }
  if(powerUnit == 'Power')
    return inpv / (Ebeam * 1.602e-10);
  else
    return inpv * 1.e-3 / 1.602e-19;
}

// Set shielding material from select
function ShieldMaterial() {
  var name = document.getElementById('NamShl').value;
  switch(name) {
    case 'CONC':
      IDShld = 0;
      break;
    case 'IRON':
      IDShld = 1;
      break;
    case 'LEAD':
      IDShld = 2;
      break;
    case 'MISC':
      IDShld = 3;
      break;
  }
}

// Input validation
function validate() {
  Ebeam = parseFloat(document.getElementById('Ebeam').value);
  if(Ebeam < 0.01) {
    Ebeam = 0.01
    document.getElementById('Ebeam').value = '0.01'
  }
  powerUnit = document.getElementById('powerMode').value;
  eps = GetEPS();

  if(document.getElementById('unit_SI').checked) {
    units = 'unit_SI';
    doseUnit = 'uSv/h';
    dosenorm = eps * 1e4 * 3600;
    lunit = cm;
    lunitname = 'cm';
  }
  else {
    units = 'unit_US';
    doseUnit = 'mrem/h';
    dosenorm = eps * 1e3 * 3600;
    lunit = inch;
    lunitname = 'inch';
  }

  IDtar = parseInt(document.getElementById('NamTar').value);
  TarLen = parseFloat(document.getElementById('TarLen').value) * lunit;
  TarRad = parseFloat(document.getElementById('TarRad').value) * lunit;

  IDShld = parseInt(document.getElementById('NamShl').value);
  AngShl = parseFloat(document.getElementById('AngShl').value);
  Theta0 = Math.ceil(AngShl - 90 + minTheta);
  Theta1 = Math.floor(AngShl + 90 - minTheta);
  Ntheta = Theta1 - Theta0;
  AngShl *= degree;
  Theta0 *= degree;
  Theta1 *= degree;

  DisShl = parseFloat(document.getElementById('DisShl').value) * lunit;
  ThkShl = parseFloat(document.getElementById('ThkShl').value) * lunit;
  AddCON = parseFloat(document.getElementById('AddCON').value) * lunit;
  AddFE  = parseFloat(document.getElementById('AddFE').value)  * lunit;
  AddPB  = parseFloat(document.getElementById('AddPB').value)  * lunit;
  AddMIS = parseFloat(document.getElementById('AddMIS').value) * lunit;

  FTDSW = parseInt(document.getElementById('FTDSW').value);
  IattSW = document.getElementById('IattSW').checked;

  return true;
}

function updateDose() {
  var tbl = document.getElementById('dose_table').tBodies[0];
  // delete all but the head row
  for (var i = tbl.rows.length - 1; i >= 0; i--) {
    tbl.deleteRow(i);
  };
  // write new dose data
  for (var i = 0; i < dose.length; i++) {
    var row = tbl.insertRow(-1)
    var cell = row.insertCell(0);
    cell.innerHTML = (Theta0/degree + i).toFixed(0)
    for (var j = 0; j < 8; j++) {
      cell = row.insertCell(j+1);
      if (dose[i][j] < 0.01) {
        cell.innerHTML = dose[i][j].toExponential(2);        
      } else {
        cell.innerHTML = dose[i][j].toPrecision(3);
      }
      
    };
  };
}

//===========================================================
// Drawing functions
//===========================================================

// Draw dose curve
function DrawDose() {
  var drawdata=[];
  for (var i = 0; i <= Ntheta; i++) {
    var Theta = Theta0 + (Theta1 - Theta0) * i / Ntheta;
    drawdata.push([Theta/degree, dose[i][5], dose[i][6], dose[i][7]]);
  };

  // Dygraph
  var g2 = new Dygraph(document.getElementById("dose_draw_div"),
              drawdata,
              {
                legend: 'always',
                // title: 'Dose',
                labels: ['Angle', 'Neutron', 'Photon', 'Total'],
                xlabel: 'Angle (degree)',
                ylabel: 'Dose (' + doseUnit + ')',
                xLabelHeight: 14,
                yLabelWidth: 14,
                // showRangeSelector: true,
                xRangePad: 5,
                axes: {
                  x: {
                    axisLabelFontSize: 12,
                    valueFormatter: function(a) { return 'Angle:' + a.toFixed(0); }
                  },
                  y: {
                    axisLabelFontSize: 12,
                    axisLabelFormatter: function(y) { return y.toPrecision(2); }
                  }
                }
              });
}


//===========================================================
// Shield11 Dose Calculation
//===========================================================
// Material parameters, in order of 'CONC',  'FE  ', 'PB  ','MISC'
var Zmat   = [13.0,  26.0,   82.0,   74.0];
var Amat   = [26.98, 55.85, 207.19, 183.85];
var RHOmat = [2.35,  7.87,  11.35,  19.30];   // density in g/cm3
var RLmat  = [26.7, 13.84,   6.37,   6.76];   // Radiation length in cm
var XMUmat = [11.1,  10.7,   14.2,   14.1];   // Moliere length in cm
var Xmfp   =[[30.0,  47.0,   97.0,   97.0],   // MFP (g/cm2) - GRN
             [55.0, 145.0,  200.0,  200.0],   //             - MID
            [120.0, 145.0,  200.0,  200.0],   //             - HEN
             [42.0,  33.6,   24.0,   25.0],   //             - GamD
            [120.0, 145.0,  200.0,  200.0]];  //             - GamI

// Data fit to HEN cross section
var E_HEN = [0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21,
             0.22, 0.23, 0.24, 0.25, 0.27, 0.30, 0.35,
             0.40, 0.45, 0.50, 0.55, 0.60, 0.70, 0.80, 0.90, 1.00];
var CS_HEN= [0.00444,0.00711, 0.0111, 0.0156, 0.0222, 0.0298, 0.0382,
             0.0489,0.0547,0.0622,0.0711,0.0889,0.116,0.162,
             0.211,0.276,0.338,0.404,0.502,0.601,0.711,0.813,1.0];

// Threshods and minimums
var ThrHEN  = 0.150      // Threshold energy (GeV) for HEN production
var TarLenG = 0.01       // Minimum target length for photons (r.l.)
var TarLenN = 17.332196  // Minimum target length for neutrons (r.l.)
var TarRadG = 1.189869   // Minimum target radius for photons
var TarRadN = 3.736411   // Minimum target radius for neutrons (Mol.units)
 
function calculate() {
  dose = [];
  for(var ia = 0; ia <= Ntheta; ia++) {  // loop for angles
    Theta = Theta0 + (Theta1 - Theta0) * ia / Ntheta;
    var AMT = AngShl - Theta;
    var CosAMT = Math.cos(AMT);
    var Rsq = Math.pow((DisShl + ThkShl) / CosAMT, 2);

    var SltShl = ThkShl * RHOmat[IDShld] / CosAMT;
    var SltCON = AddCON * RHOmat[0] / CosAMT;
    var SltFE  = AddFE  * RHOmat[1] / CosAMT;
    var SltPB  = AddPB  * RHOmat[2] / CosAMT;
    var SltMIS = AddMIS * RHOmat[3] / CosAMT;

    var da = [];
    for(var isrc = 0; isrc < 5; isrc++)
      da.push( Src(isrc, Theta) / Rsq *
           Math.exp(-SltShl / Xmfp[isrc][IDShld]) *
           Math.exp(-SltCON / Xmfp[isrc][0]) *
           Math.exp(-SltFE  / Xmfp[isrc][1]) *
           Math.exp(-SltPB  / Xmfp[isrc][2]) *
           Math.exp(-SltMIS / Xmfp[isrc][3]) * dosenorm);
    da.push(da[0] + da[1] + da[2]);  // [5] Neutron dose
    da.push(da[3] + da[4]);          // [6] Photon dose
    da.push(da[5] + da[6]);          // [7] Total dose
    da.push(da[6] / da[5]);          // [8] Photon/Neutron ratio
    dose.push(da);
  }
}

function Src(type, angle) {
  var AbsAng = Math.abs(angle);
  var CosThe = Math.cos(AbsAng);
  var SinThe = Math.sin(AbsAng);
  var TarLenGCM = TarLen * RHOmat[IDtar];
  var TarLenRL  = TarLenGCM / RLmat[IDtar];
  var TarRadGCM = TarRad * RHOmat[IDtar];
  var TarRadMU  = TarRadGCM / XMUmat[IDtar];
  var TarRadRelax = TarRadGCM / Xmfp[3][IDtar];
  var SltSor = 0.0;
  var CritTar, SltTar, CritCor, SltCor, Sorc;

  // Check for neutron attenuation by target (GRN, MID, HEN, and GamI only)
  if(IattSW && type != 4) {
    if(AbsAng > Math.PI/2) {
      SltSor = 0.0;
    }
    else {
      // Find critical angle and SltTar for cylindrical target
      CritTar = Math.atan(TarRad / TarLen)
      if(AbsAng >= CritTar)
        SltTar = TarRadGCM / SinThe;
      else
        SltTar = TarLenGCM / CosThe;
      // Find critical angle and SltCor for cylindrical core
      CritCor = Math.atan(TarRadN * XMUmat[IDtar] / TarLenN / RLmat[IDtar]);
      if(AbsAng >= CritCor)
        SltCor = TarRadN * XMUmat[IDtar] / SinThe;
      else
        SltCor = TarLenN * RLmat[IDtar] / CosThe;

      SltSor = SltTar - SltCor;
    }
  }
  switch(type) {
    case 0:  // GRN: Giant-Resonance Neutrons
      if(FTDSW == 0)
        Sorc = 4.93 * Math.pow(Zmat[IDtar], 0.662);
      else if(FTDSW == 1)
        Sorc = 4.93 * Math.pow(Zmat[IDtar], 0.662) * (3.6e-8/3.2e-8);
      else if(FTDSW == 2)
        Sorc = 4.93 * Math.pow(Zmat[IDtar], 0.662) * (4.69e-8/3.2e-8);
      break;
    
    case 1:  // MID: Mid-Energy Neutrons
      if(FTDSW == 0)
        Sorc = 43.9 / Math.pow(Amat[IDtar], 0.37) / (1.0 - 0.75 * CosThe);
      else if(FTDSW == 1)
        Sorc = 43.9 / Math.pow(Amat[IDtar], 0.37) / (1.0 - 0.75 * CosThe) * (4.5e-8/3.2e-8);
      else if(FTDSW == 2)
        Sorc = 43.9 / Math.pow(Amat[IDtar], 0.37) / (1.0 - 0.75 * CosThe) * (9.94e-8/3.2e-8);
      
      if(Ebeam <= 0.5)
        Sorc *= 1.6 * Math.pow(Ebeam, 1.5);
      else if(Ebeam > 0.5 && Ebeam < 1.0)
        Sorc *= 0.566 + 0.434 * (Ebeam - 0.5) / 0.5;
      break;

    case 2:  // HEN
      if(Ebeam <= ThrHEN)
        Sorc = 0.0
      else {
        if(FTDSW == 0)
          Sorc = 13.7 / Math.pow(Amat[IDtar], 0.65) / Math.pow((1.0 - 0.72 * CosThe), 2);
        else if(FTDSW == 1)
          Sorc = 13.7 / Math.pow(Amat[IDtar], 0.65) / Math.pow((1.0 - 0.72 * CosThe), 2) * (6.5e-8/6.7e-8);
        else if(FTDSW == 2)
          Sorc = 13.7 / Math.pow(Amat[IDtar], 0.65) / Math.pow((1.0 - 0.72 * CosThe), 2) * (1.23e-8/6.7e-8);
        if(Ebeam < 1.0) {
          for(var k = 1; k < 23; k++) {
            if(Ebeam < E_HEN[k]) {
              var DelCS = (CS_HEN[k] - CS_HEN[k-1]) * (Ebeam - E_HEN[k-1]) / (E_HEN[k] - E_HEN[k-1]);
              Sorc *= CS_HEN[k-1] + DelCS;
              break;
            }
          }
        }
      }
      break;

    case 3:  // GammaD
      Sorc = 1.06e6 * Ebeam * Math.exp(-TarLenGCM / Xmfp[3][IDtar]) * Math.exp(-0.959 * Math.sqrt(AbsAng/degree));
      if(AbsAng <= Math.PI/2)
        Sorc += 683.0 * Math.exp(-TarRadRelax) * Math.exp(-AbsAng/degree/72.2);
      else
        Sorc += 683.0 * Math.exp(-TarRadG) * Math.exp(-AbsAng/degree/72.2);
      break;

    case 4:  // GammaI
      if(Ebeam <= ThrHEN)
        Sorc = 0.0
      else {
        if(FTDSW == 0)
          Sorc = 13.7 / Math.pow(Amat[IDtar], 0.65) / Math.pow((1.0 - 0.72 * CosThe), 2);
        else if(FTDSW == 1)
          Sorc = 13.7 / Math.pow(Amat[IDtar], 0.65) / Math.pow((1.0 - 0.72 * CosThe), 2) * (6.5e-8/6.7e-8);
        else if(FTDSW == 2)
          Sorc = 13.7 / Math.pow(Amat[IDtar], 0.65) / Math.pow((1.0 - 0.72 * CosThe), 2) * (1.23e-8/6.7e-8);
        if(Ebeam < 1.0) {
          for(var k = 1; k < 23; k++) {
            if(Ebeam < E_HEN[k]) {
              var DelCS = (CS_HEN[k] - CS_HEN[k-1]) * (Ebeam - E_HEN[k-1]) / (E_HEN[k] - E_HEN[k-1]);
              Sorc *= CS_HEN[k-1] + DelCS;
              break;
            }
          }
        }
      }
      Sorc *= 0.267;
      break;
  }
  Sorc *= Ebeam * 1.e-11;
  if(IattSW)
    Sorc *= Math.exp(-SltSor / Xmfp[type][IDtar]);
  return Sorc;
}

