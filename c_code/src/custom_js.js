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
  }
});
