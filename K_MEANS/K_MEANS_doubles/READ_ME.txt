********************************************************************************
K-means algorithm compilation: make
K-means algorithm run : ./seq_main -o -t 0.0005 -b -n 1000 -i Image_data/color17695.bin  (PRINT->OFF)
                      or  ./seq_main -o -t 0.0005 -b -d -n 1000 -i Image_data/color17695.bin  (PRINT->ON)
Check_approximation_accuracy compilation :  gcc -Wall -g check_approximation_accuracy.c -o check_approximation_accuracy -lm
Check_approximation_accuracy run : ./check_approximation_accuracy -n 1000 (PRINT->OFF)
                                or ./check_approximation_accuracy -d -n 1000 (PRINT->ON)

                                
Object sampling: 1)Change line 156 in file seq_main.c

Dimension sampling 1)Change line 121 in file seq_main.c
                   2)Set dim variable in line 148 in file seq_main.c
                   3)Change line 301 in file check_approximation_accuracy.c
