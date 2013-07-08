define(['require', 'kievII', 'image'], function(require, K2) {

  
    var imgResources = null;

    /* This gets returned to the host as soon as the plugin is loaded */ 
    var pluginConf = {
        name: "Voron",
        osc: false,
        audioIn: 1,
        audioOut: 1,
        ui: {
            type:'canvas',
            width: 268,
            height: 340
        },
    };

    /* This gets called when all the resources are loaded */
    var pluginFunction = function (args, resources) {
      
        var backgroundImage = resources[2];
        var knobImage = resources[3];
		var discLeft = resources[4];
		var discRight = resources[5];
        var nSamples = 2048;
        var fftFrameSize = 2048;
        shifterStartValue = 0;
        
        console.log ("plugin inited, args is", args, "KievII object is ", K2, "bg image is", backgroundImage);
        
        this.name = args.name;
        this.id = args.id;
        
        // The sound part
        this.audioSource = args.audioSources[0];
        this.audioDestination = args.audioDestinations[0];
        this.context = args.audioContext;
        var context = this.context;
        
        this.processorNode = this.context.createJavaScriptNode(nSamples, 1, 1);
        
        this.shifter = new Pitchshift(fftFrameSize, this.context.sampleRate, 'FFT');
        
        this.processorNode.onaudioprocess = function (event) {
            // Get left/right input and output arrays
            var outputArray = [];
            outputArray[0] = event.outputBuffer.getChannelData(0);
            var inputArray = [];
            inputArray[0] = event.inputBuffer.getChannelData(0);

            var data = inputArray[0];
            this.shifter.process (this.shiftValue, data.length, 4, data);
            
            var out_data = outputArray[0];
            for (var i = 0; i < out_data.length; ++i) {
                out_data[i] = this.shifter.outdata[i];
            }
            
        }.bind(this);
        
        this.audioSource.connect (this.processorNode);
        this.processorNode.connect (this.audioDestination);
        
        // The graphical part
        this.ui = new K2.UI ({type: 'CANVAS2D', target: args.canvas});
        
        this.viewWidth = args.canvas.width;
        this.viewHeight = args.canvas.height;
        
        /* BACKGROUND INIT */
        
        var bg = new K2.Background({
            ID: 'background',
            image: backgroundImage,
            top: 0,
            left: 0
        });
    
        this.ui.addElement(bg, {zIndex: 0});
        
        /* KNOB INIT */
       var knobArgs = {
            ID: "pitch_knob",
            left: 87 ,
            top: 176,
            image : knobImage,
            sensitivity : 5000,
            initAngValue: 270,
            startAngValue: 218,
            stopAngValue: 501,
            /* knobMethod: 'updown', */
            onValueSet: function (slot, value) {
                var shift_value = value * (1.5) + 0.5;
                /* shift argument is like a play rate */
                /* We want 0 -> 0.5, 0.5 -> 1, 1 -> 2 */
                /* Let's calculate the semitones */
                var semitoneShift =  K2.MathUtils.linearRange (0, 1, -12, 12, value);
				if (this.discrete === 1) {
					semitoneShift = Math.round(semitoneShift);
				}
                /* Let's calculate the "play rate" */
                shift_value = Math.pow(1.0595, semitoneShift);
                this.shiftValue = shift_value;
                console.log ('Shift value set to ', value, this.shiftValue);
                this.ui.refresh();
            }.bind(this),
            isListening: true
        };
        
        this.ui.addElement(new K2.RotKnob(knobArgs));
        this.ui.setValue({elementID: "pitch_knob", value: 0.5});

		/* Button init */
        var buttonArgs = {
            ID: "discButton",
            left: 104,
            top: 108,
            imagesArray : [discLeft, discRight],
            onValueSet: function (slot, value) {
				this.discrete = value;
                this.ui.refresh();
            }.bind(this),
            isListening: true
        };
        
        this.ui.addElement(new K2.Button(buttonArgs));

        this.ui.refresh();

        // Initialization made it so far: plugin is ready.
        args.hostInterface.setInstanceStatus ('ready');
    };
  
    /* This function gets called by the host every time an instance of
       the plugin is requested [e.g: displayed on screen] */        
    var initPlugin = function (initArgs) {
        var args = initArgs;

        var requireErr = function (err) {
            var failedId = err.requireModules && err.requireModules[0];
            require.undef(failedId);
            args.hostInterface.setInstanceStatus ('fatal', {description: 'Error initializing plugin: ' + failedId});
        }.bind(this);

        var resList = [ 'https://github.com/corbanbrook/dsp.js/raw/master/dsp.js',
                        'https://github.com/janesconference/KievII/raw/master/dsp/pitchshift.js',
                        'image!./assets/images/Voron_bg2.png!rel',
                        'image!./assets/images/white_big.png!rel',
                        'image!./assets/images/switch_l.png!rel',
                        'image!./assets/images/switch_r.png!rel'];
        
        require (resList,
            function () {
                var resources = arguments;
                pluginFunction.call (this, args, resources);
            }.bind(this),
            function (err) {
                requireErr (err);
            });
    };
    
    return {
        initPlugin: initPlugin,
        pluginConf: pluginConf
    };
});