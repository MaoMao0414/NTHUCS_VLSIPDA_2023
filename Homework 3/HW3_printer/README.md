# Homework 3 Printer
You can use this program to visualize your floorplanning result.

## Preparing
* **Step1:**  
    This floorplan visualization tool is written in **Processing**. Please download the [**Processing**](https://processing.org/download) and unzip it.

* **Step2:**  
    Unzip `HW3_printer.zip` and put this directory in the folder `processing`. The file hierarchy is shown below:
    ```
    processing/
    |   core/
    |   HW3_printer/
    |   |   figure/
    |   |   output/
    |   |   testcase/
    |   |   |   public1.txt
    |   |   |   public2.txt
    |   |   |   ...
    |   |   HW3_printer.pde
    |   |   README.md
    |   java/
    |   lib/
    |   modes/
    |   tools/
    |   processing-java.exe
    |   processing.exe
    |   revisions.md
    ```

* **Step3:**  
    Put your floorplanning result in the folder `output`.
    ```
    processing/
    |   core/
    |   HW3_printer/
    |   |   figure/
    |   |   output/
    |   |   |   public1.floorplan
    |   |   |   public2.floorplan
    |   |   |   ...
    ```

## Working Flow
* **Step1:**  
    Open `HW3_printer.pde` with `processing.exe`. You can change the testcase by changing the *caseName* in **setup()**.

    ```java
    void setup() {
        size(800, 800);              // image size
        String caseName = "public1"; // change testcase

        Printer printer = new Printer();
        printer.readTxt("./testcase/" + caseName + ".txt");
        printer.readFloorplan("./output/" + caseName + ".floorplan");
        printer.drawFloorplan("./figure/" + caseName + ".png");
    }
    ```
    The default image size is 800 x 800. If you want to change the image size, please revise the value of **size(x, x)**.

* **Step2:**  
    Press the **run** button to execute the drawing program. A drawing window will pop up on your screen. The floorplanning image will be saved in the folder `figure`.
    ```
    processing/
    |   HW3_printer/
    |   |   figure/
    |   |   |   public1.png
    |   |   |   public2.png
    |   |   |   ...
    ```

## Reference
* https://processing.org/reference/
