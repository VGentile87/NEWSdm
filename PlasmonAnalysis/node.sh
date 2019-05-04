#!bin/bash
root -l debug6_test_grain.root <<EOC
tree1->MakeClass("myNode")
.q
EOC
rm myNode.C
