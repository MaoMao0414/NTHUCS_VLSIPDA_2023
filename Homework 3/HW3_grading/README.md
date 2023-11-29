# Homework 3 Grading Script
We will use this script to grade your program. **Make sure your program can be executed by this script.**

## Preparing
* Step1:  
    Enter the `HW3_grading` directory and create a new directory named with your student ID in the `student` directory.
    ```sh
    $ cd HW3_grading/
    $ mkdir student/${your_student_id}
    ```
    For example,
    ```sh
    $ cd HW3_grading/
    $ mkdir student/112062500
    ```

* Step2:  
    Put your compressed file in the directory which you just created.  
    The correct path should be:
    ```
    HW3_grading/student/${your_student_id}/CS6135_HW3_${your_student_id}.tar.gz
    ```
    For example,
    ```
    HW3_grading/student/112062500/CS6135_HW3_112062500.tar.gz
    ```

### Notice:
**Please make sure not to put your original directory here**, as it will remove all directories before unzipping the compressed file.

## Grading
* Step1:  
    Navigate to the `HW3_grading` directory and run `HW3_grading.sh`.
    ```sh
    $ cd HW3_grading/
    $ bash HW3_grading.sh
    ```

* Step2:  
    Check your output. Ensure the status of each checking item is **yes**.
    * If the status of a testcase is **success**, it means your program finished in time, and the output result is legal.
        ```
        grading on 112062500:
         checking item          | status
        ------------------------|--------
         correct tar.gz         | yes
         correct file structure | yes
         have README            | yes
         have Makefile          | yes
         correct make clean     | yes
         correct make           | yes

          testcase | wirelength |    runtime | status
        -----------|------------|------------|--------
           public1 |  157924052 |     590.63 | success
           public2 |   19610550 |     591.21 | success
           public3 |    1882293 |     591.42 | success
           public4 |   63614225 |     591.80 | success
        ```

    * If the status of a testcase is not **success**, it means your program failed in this testcase.
        ```
        grading on 112062500:
         checking item          | status
        ------------------------|--------
         correct tar.gz         | yes
         correct file structure | yes
         have README            | yes
         have Makefile          | yes
         correct make clean     | yes
         correct make           | yes

          testcase | wirelength |    runtime | status
        -----------|------------|------------|--------
           public1 |        N/A |        TLE | Time out while testing sample.
           public2 |        N/A |        TLE | Time out while testing sample.
           public3 |        N/A |        TLE | Time out while testing sample.
           public4 |        N/A |        TLE | Time out while testing sample.
        ```
