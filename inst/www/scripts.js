/*
	Toggles disabled property of the generateBtn

	param disabled Boolean
*/
function toggleGenerateBtn(disabled) {
	$('#generateBtn').prop('disabled', disabled);
}

/*
	Toggle the visibility of the progress overlay

	param running Boolean
*/
function toggleProgress(running) {
	$('#progressOverlay').toggleClass('hidden', !running);
}

/*
	Toggles the showing of the success state of generation

	param success Boolean
*/
function toggleSuccess(success) {
	$('#progressOverlay .progress').toggleClass('active', !success);
	$('#progressOverlay .bar').toggleClass('bar-success', success);
	$('#progressOverlay a.success').toggleClass('hidden', !success);
}

/*
	Toggles the showing of error state of generation

	param error Boolean
	param message String
*/
function toggleError(error) {
	$('#progressOverlay .progress').toggleClass('active', !error.state);
	$('#progressOverlay .bar').toggleClass('bar-danger', error.state);
	$('#progressOverlay button.failure').toggleClass('hidden', !error.state);
	$('#progressOverlay .errorMessage')
		.text(error.message)
		.toggleClass('hidden', !error.state);
}

/*
	Removes progress overlay and reset file upload
*/
function retryUpload() {
	Shiny.unbindAll();

	toggleProgress(false);
	toggleSuccess(false);
	toggleError({state: false, message: ''});
	$('#panSelect').replaceWith($('#panSelect').get(0).outerHTML);
	$('#panSelect').on('change', function() {toggleGenerateBtn(false)});
	toggleGenerateBtn(true);

	Shiny.bindAll();

}

$(document).ready(function() {
	// Shiny handlers
	Shiny.addCustomMessageHandler('toggleGenerateBtn', toggleGenerateBtn);
	Shiny.addCustomMessageHandler('toggleProgress', toggleProgress);
	Shiny.addCustomMessageHandler('toggleSuccess', toggleSuccess);
	Shiny.addCustomMessageHandler('toggleError', toggleError);

	// UI handlers
	$('#progressOverlay button.failure').on('click', retryUpload);
	$('#panSelect').on('change', function() {toggleGenerateBtn(false)});
	$('#generateBtn').on('click', function() {toggleProgress(true)});
	$('#download').on('click', retryUpload);
})