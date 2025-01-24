awk 'BEGIN {atom_count=1} 
     /^ATOM|^HETATM/ {
         printf "%s%5d%s\n", substr($0, 1, 6), atom_count, substr($0, 12);
         atom_count++;
     }
     !/^ATOM|^HETATM/ { print }'  $1
