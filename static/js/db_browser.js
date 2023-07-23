        $(document).ready(function() {
            $("#edit-button").click(function() {
                $(".editable-field").prop("disabled", false);
                $(".submit-button").prop("disabled", false);
                $(this).prop("disabled", true);
            });
        });

function goToHomepage() {
    window.location.href = "/";
    }