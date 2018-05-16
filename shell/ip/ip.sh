#/bin/bash
#
# The purpose of the script is to see which IP is occupied in the current LAN.

a=0
echo >ip.txt
while :
do
    a=$(($a+1))
    if test $a -gt 255
    then break
    else
        echo $(ping -c 1 192.168.1.$a | grep "ttl" | awk '{print $4}'| sed 's/://g')
        ip=$(ping -c 1 192.168.1.$a | grep "ttl" | awk '{print $4}'| sed 's/://g')
        echo $ip >> ip.txt &
    fi
done
sed -i '/^$/d' ip.txt
