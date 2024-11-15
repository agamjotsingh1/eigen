#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "Incorrect usage of shell script. You must pass exactly one argument, i.e., the number of test cases."
    exit 1
fi

input_file=./tests/test.c
verify_file=./tests/verify.py
num_testcases=$1

# Check if gcc is installed
if command -v gcc &> /dev/null
then
    echo "gcc is installed"
else
    echo "Install GCC"
    echo "Exiting..."
    exit 1
fi

# Compile the input file
echo "Compiling..."
gcc $input_file -o $input_file.out -lm
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi
echo "Compilation successful!"

test(){
    testNo=$1
    inFile="tests/input/input_$testNo.txt"
    outFile="tests/output/output_$testNo.txt"
    expectedFile="tests/expected/expected_$testNo.txt"

    cat $inFile | "$input_file.out" > $outFile
    cat $inFile | python3 "$verify_file" > $expectedFile
    res=`diff $outFile $expectedFile | wc -l`
    if [ $res -eq 0 ]
    then
        echo "TESTCASE-$testNo PASSED"
    else
        echo "TESTCASE-$testNo FAILED"
    fi
}

for i in $(seq 1 $num_testcases)
do
    test $i
done
