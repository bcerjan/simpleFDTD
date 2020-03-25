// Custom JS library to merge in functions we need to call from C
mergeInto(LibraryManager.library, {
  updateProgress: function(percent) {
    percent = percent.toFixed(3);
    const bar = document.querySelector('#progress');
    bar.innerHTML = percent+"%";
    bar.style.width = percent + "%";
    bar.setAttribute('aria-valuenow', percent);
    console.log('Inside JS Func!');
  }
});
