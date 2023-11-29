/* ||=========================================|| */
/* || PDA HW3 floorplan visualization program || */
/* ||=========================================|| */

void setup() {
    size(800, 800);              // image size
    String caseName = "public1"; // change testcase

    Printer printer = new Printer();
    printer.readTxt("./testcase/" + caseName + ".txt");
    printer.readFloorplan("./output/" + caseName + ".floorplan");
    printer.drawFloorplan("./figure/" + caseName + ".png");
}

class Module {
    String name;
    float x = 0, y = 0, w = 0, h = 0;
    boolean isSoft = false;
    
    Module(String name, float x, float y, float w, float h, boolean isSoft) {
        this.name = name;
        this.x = x;
        this.y = y;
        this.w = w;
        this.h = h;
        this.isSoft = isSoft;
    }
}

class Printer {
    float chipWidth = 0, chipHeight = 0;
    ArrayList<Module> modules = new ArrayList<Module>();
    float scaleFactor = 1;

    void readTxt(String filename) {
        BufferedReader fin = createReader(filename);
        String line = null;
        try {
            while ((line = fin.readLine()) != null) {
                String[] tokens = splitTokens(line);
                if (tokens.length == 0)
                    continue;
                
                if (tokens[0].equals("ChipSize")) {
                    chipWidth = float(tokens[1]);
                    chipHeight = float(tokens[2]);
                }
                else if (tokens[0].equals("FixedModule")) {
                    modules.add(new Module(tokens[1], float(tokens[2]), float(tokens[3]),
                        float(tokens[4]), float(tokens[5]), false));
                }
            }
            fin.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    void readFloorplan(String filename) {
        BufferedReader fin = createReader(filename);
        String line = null;
        try {
            while ((line = fin.readLine()) != null) {
                String[] tokens = splitTokens(line);
                if (tokens.length == 0)
                    continue;
                
                if (tokens.length == 5) {
                    modules.add(new Module(tokens[0], float(tokens[1]), float(tokens[2]),
                        float(tokens[3]), float(tokens[4]), true));
                }
            }
            fin.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    void drawFloorplan(String filename) {
        float scaleByWidth = (width - 1) / chipWidth;
        float scaleByHeight = (height - 1) / chipHeight;
        scaleFactor = (scaleByWidth < scaleByHeight)? scaleByWidth : scaleByHeight;

        background(255);
        drawModules();
        drawOutline();
        save(filename);
    }
    
    float scaleVal(float val) {
        return val * scaleFactor;
    }
    
    // vertical flip
    float getY(Module module) {
        return chipHeight - (module.y + module.h);
    }
    
    void drawModules() {
        for (Module module : modules) {
            float x = scaleVal(module.x);
            float y = scaleVal(getY(module));
            float w = scaleVal(module.w);
            float h = scaleVal(module.h);
            stroke(0);
            if (module.isSoft)
                fill(200);
            else
                fill(100);
            rect(x, y, w, h);
            
            fill(0);
            textSize(12);
            text(module.name, x + 2, y + 14);
        }
    }

    void drawOutline() {
        float w = scaleVal(chipWidth);
        float h = scaleVal(chipHeight);
        stroke(255, 0, 0);
        strokeWeight(1.5);
        noFill();
        rect(0, 0, w, h);
    }
}
