
export DYLD_LIBRARY_PATH=$HOME/local/lib

while read -r line
do
    echo "====================="
    echo "Running $line"
    echo "====================="
    
    ./fact "../../Netlib/data/$line" 1 1
done < "../../Netlib/netlib_names.txt"
