if [ $# -eq 2 ]
then args="$1 $2"
else args=$1
fi
echo $args

make && time src/meshpro $args forest1.off+ && src/meshview forest1.off+
	
