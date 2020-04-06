/**
    Copyright (c) 2020 Ben Cerjan

    This file is part of simpleFDTD.

    simpleFDTD is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    simpleFDTD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with simpleFDTD.  If not, see <https://www.gnu.org/licenses/>.
**/

// Custom JS library to merge in functions we need to call from C

/*
   Note that this file has hard-coded id's for various things it references
   so if you change simulation_template.html you might also need to adjust
   the naming of things in this file.
*/

// Assumes progess bar has id="progress"
mergeInto(LibraryManager.library, {
  updateProgress: function(percent) {
    percent = percent.toFixed(3);
    const bar = document.querySelector('#progress');
    bar.innerHTML = percent+"%";
    bar.style.width = percent + "%";
    bar.setAttribute('aria-valuenow', percent);
  }
});


// See: https://www.chartjs.org/docs/latest/developers/updates.html
// Assumes our data chart has var outputChart = new Chart({...})
// and that refl data is the first dataset and tran data is the second
mergeInto(LibraryManager.library, {
  addReflData: function(x,y) {
    data = {
      x: x,
      y: y
    };
    outputChart.data.datasets[0].data.push(data);
  }
});

mergeInto(LibraryManager.library, {
  addTranData: function(x,y) {
    data = {
      x: x,
      y: y
    };
    outputChart.data.datasets[1].data.push(data);
  }
});

// Split this into two functions so we add data and then only update the chart once
// Assumes our data chart has var outputChart = new Chart({...})
mergeInto(LibraryManager.library, {
  updateChartData: function() {
    outputChart.update();
    console.log(outputChart);
  }
});

// Function to re-enable "Run Simulation" button after we finish:
mergeInto(LibraryManager.library, {
  enableButton: function() {
    runSim.disabled = false;
    runSim.innerHTML = "Run Simulation";
  }
});

// Function to re-enable "Run Simulation" button after we finish:
mergeInto(LibraryManager.library, {
  enableButton: function() {
    runSim.disabled = false;
    normData.disabled = false;
    runSim.innerHTML = "Run Simulation";
  }
});
