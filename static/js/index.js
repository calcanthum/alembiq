        function runDBInit() {
            // Sending a GET request to the /db_init endpoint
            fetch('/db_init')
                .then(response => response.text())  // Change to text() to handle HTML response
                .then(data => {
                    console.log(data);
                    // Display the response in a dialog box
                    alert(data);
                })
                .catch(error => console.error('Error:', error));
        }

        function confirmDBInit() {
            // Display a confirmation dialog to the user
            if (confirm("Are you sure you want to initialize the database? This action is irreversible.")) {
                // If the user clicks "OK," proceed with the database initialization
                runDBInit();
            } else {
                // If the user clicks "Cancel," do nothing
            }
        }

        function runMolDraw() {
            // Redirecting to the /moldraw endpoint
            window.location.href = '/moldraw';
        }

        function insertMolecule() {
            const molecule_string = prompt('Please enter the SMILES or InChI string:');
            if (molecule_string) {
                // Sending a POST request to the /molecule_insert endpoint with the molecule_string
                fetch('/molecule_insert', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json'
                    },
                    body: JSON.stringify({ molecule_string: molecule_string })
                })
                .then(response => response.json())
                .then(data => {
                    console.log(data);
                    // Display the response in a dialog box
                    alert(data.message);
                })
                .catch(error => console.error('Error:', error));
            }
        }

        function runDBBrowser() {
            // Redirecting to the /db_browser endpoint
            window.location.href = '/db_browser';
        }