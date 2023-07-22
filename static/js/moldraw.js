var viewer = null;

function displayStructure(data) {
    if (viewer === null) {
        viewer = new NGL.Stage('molecule-viewer');
    } else {
        viewer.removeAllComponents();
    }
    viewer.loadFile(new Blob([data.pdb_data], {type: 'text/plain'}), { ext: 'pdb' })
    .then(function(comp) {
        changeRepresentation();
        viewer.autoView();
    });
}

function fetchMolecule(url, options = {}) {
    fetch(url, options)
    .then(response => response.json())
    .then(data => displayStructure(data))
    .catch(error => console.error('Error:', error));
}

function displayInputStructure() {
    var moleculeString = document.getElementById('input-field').value;
    fetchMolecule('/molecule', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ molecule_string: moleculeString })
    });
    document.getElementById('input-field').value = '';
}

function computeProperties() {
    var moleculeString = document.getElementById('input-field').value;
    fetch('/compute_properties', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({molecule_string: moleculeString}),
    })
    .then(response => response.json())
    .then(data => {
        var propertiesText = `Molecular Weight: ${data.molecular_weight}\n` +
            `Number of H Donors: ${data.num_h_donors}\n` +
            `Number of H Acceptors: ${data.num_h_acceptors}`;
        document.getElementById('computed-properties').value = propertiesText;
    })
    .catch((error) => {
        console.error('Error:', error);
    });
}

function changeRepresentation() {
    var representation = document.querySelector('input[name="representation"]:checked').value;
    if (viewer) {
        viewer.eachComponent(function(comp){
            comp.removeAllRepresentations();
            comp.addRepresentation(representation);
        });
    }
}

window.addEventListener('load', function() {
    fetch('/molecule_identifiers')
    .then(response => response.json())
    .then(identifiers => {
        var selector = document.getElementById('molecule-selector');
        identifiers.forEach(([id, smiles]) => {
            var option = document.createElement('option');
            option.value = id;
            option.text = smiles;
            selector.appendChild(option);
        });

        // Add event listeners to all radio buttons that change the representation
        document.querySelectorAll('input[name="representation"]').forEach(function(radio) {
            radio.addEventListener('change', function() {
                if (viewer) {
                    changeRepresentation(); // Call changeRepresentation() when the radio button selection changes
                }
            });
        });
    })
    .catch(error => console.error('Error:', error));
});


function displaySelectedStructure() {
    var id = document.getElementById('molecule-selector').value;
    fetchMolecule('/molecule/' + id);
}

// Use the 'change' event instead of 'click' for the radio buttons
document.querySelectorAll('input[name="representation"]').forEach(function(radio) {
    radio.addEventListener('change', function() {
        changeRepresentation(); // Call changeRepresentation() when the radio button selection changes
    });
});