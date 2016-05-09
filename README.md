# projetcpp

https://liris.cnrs.fr/~crocaber/3Bim/Projet/

set pointsize 1
plot "< awk '{if($3 == \"0\") print}' output1_1500_0.9.txt" u 1:2 t "red" w p pt 2, \
     "< awk '{if($3 == \"1\") print}' output1_1500_0.9.txt" u 1:2 t "green" w p pt 2, \
     "< awk '{if($3 == \"2\") print}' output1_1500_0.9.txt" u 1:2 t "blue" w p pt 2, \
     "< awk '{if($3 == \"3\") print}' output1_1500_0.9.txt" u 1:2 t "orange" w p pt 2

set xrange [1:1500]
set term png
set output 'd001.png'
set palette model RGB defined (0 "salmon",1 "web-blue", 2 "light-green")
plot 'd0001.txt' using 1:2:3 notitle with points pt 157 palette
