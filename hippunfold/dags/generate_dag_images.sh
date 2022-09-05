#pipe hippunfold into this command


if [ "$#" -lt 1 ]
then
    echo "Pipe hippunfold with --dag or --rulegraph into this command"
    echo "Usage: $0 <name> "
    exit 1
fi

name=$1

cat | ./add_subgraphs.sh > txt/${name}.txt
cat txt/${name}.txt | dot -Tpdf > pdf/${name}.pdf
cat txt/${name}.txt | dot -Tsvg > svg/${name}.svg 


