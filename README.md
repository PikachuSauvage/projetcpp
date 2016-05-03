# projetcpp

https://liris.cnrs.fr/~crocaber/3Bim/Projet/

set pointsize 1
plot "< awk '{if($3 == \"0\") print}' output.txt" u 1:2 t "red" w p pt 2, \
     "< awk '{if($3 == \"1\") print}' output.txt" u 1:2 t "green" w p pt 2, \
     "< awk '{if($3 == \"2\") print}' output.txt" u 1:2 t "blue" w p pt 2, \
     "< awk '{if($3 == \"3\") print}' output.txt" u 1:2 t "orange" w p pt 2
