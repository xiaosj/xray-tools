// Constant
var e  = 1.602E-19;  // electron charge
var m0 = 9.109E-31;  // electron mass
var c  = 2.998E+08;  // speed of light
var h  = 6.626E-34;  // Plank constant
var pi = 3.14159;
var out_digits = 4;  // number of digits in output

// facilitate functions
function getEnergy() {   // Electron energy in GeV
  return parseFloat(document.getElementById('e_energy').value);
}

function getPeriod() {   // Undulator period in mm
  return parseFloat(document.getElementById('period_mm').value);
}

function getK() {
  return parseFloat(document.getElementById('K').value);
}

function getUndB() {
  return parseFloat(document.getElementById('und_B').value);
}

// Functions for bending magnet
function BendCriticalEnergy(Energy_GeV, B) {
  return 0.665 * Energy_GeV * Energy_GeV * B;
}

// Functions for undulaotrs
function K2B(K, period_mm) {  // Calculate B from K
  return K * m0 * c * 2 * pi / (e * period_mm / 1000);
}

function B2K(B, period_mm) {  // Calculate K from B
  return B * e * period_mm / 1000 / (m0 * c * 2 * pi);
}

function gammaE(Energy_GeV) {  // Lorentz factor of electron
  return Energy_GeV * 1000 / 0.511;
}

function und_fundamental_wavelength_nm(Energy_GeV, period_mm, K) {
  var gamma = gammaE(Energy_GeV);
  return period_mm / (2 * gamma * gamma) * (1 + K * K / 2) * 1e6;
}

function und_fundamental_energy_keV(Energy_GeV, period_mm, K) {
  var wl = und_fundamental_wavelength_nm(Energy_GeV, period_mm, K);
  return h * c / (wl * 1.e-9) / e / 1000;
}

function und_total_sr_power_mA_m(Energy_GeV, B) {
  return 0.6328 * Energy_GeV * Energy_GeV * B * B;
}

// Functions to update values
function updateUndB() {
  document.getElementById('und_B').value = K2B(getK(), getPeriod()).toFixed(out_digits);
}

function updateUndK() {
  // var B = parseFloat(document.getElementById('und_B').value);
  document.getElementById('K').value = B2K(getUndB(), getPeriod()).toFixed(out_digits);
}

function updateUndSR() {
  var energy = getEnergy();
  var period = getPeriod();
  var K = getK();
  document.getElementById('fundamental_wavelength').value = und_fundamental_wavelength_nm(energy, period, K).toFixed(out_digits);
  document.getElementById('fundamental_energy').value = und_fundamental_energy_keV(energy, period, K).toFixed(out_digits);
  document.getElementById('und_total_sr_power').value = und_total_sr_power_mA_m(energy, getUndB()).toFixed(out_digits);
}

function updateBendSR() {
  var energy = getEnergy();
  var B = parseFloat(document.getElementById('bend_B').value);
  document.getElementById('bend_CE').value = BendCriticalEnergy(energy, B).toFixed(out_digits);
}

function init() {
  updateBendSR();
  updateUndB();
  updateUndSR();
}

function changeEnergy() {
  updateBendSR();
  updateUndSR();
}

function changeBendB() {
  updateBendSR();
}

function changeUndB() {
  updateUndK();
  updateUndSR();
}

function changeK() {
  updateUndB();
  updateUndSR();
}

function changePeriod() {
  updateUndB();
  updateUndSR();
}
