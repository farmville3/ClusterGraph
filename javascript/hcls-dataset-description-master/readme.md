HCLS Dataset Descriptions statistics scripts
====

This is a work in progress.

## Statistics - hcls-dd.rb

Scripts to produce statistics over a SPARQL endpoint that are compliant with the W3C Dataset Descriptions for Health Care and Life Sciences (HCLS) as provided by the HCLS interest group. See https://www.w3.org/TR/hcls-dataset/

Note that the script will compute the statistics for each graph in the endpoint individually.

```
ruby hcls-dd.rb


    hcls-dd.rb [options] <SPARQL url>

    Input:
       --uri               <endpoint uri>
       --sparql            <target sparql uri>
       --default-graph     <default graph uri>
       --label             <,'endpoint label'> [default: my_endpoint]
       --description       <'endpoint description'> [default: My Endpoint Statistics]

    Output:
       --output_json       <outfile.json> save hash as json file
       --output_rdf        <outfile.ttl>

    Options:
       --basic_metrics     to compute only the basic metrics

       --help              print this

```

The script will produce a ntriple graph file (for the dataset description) and a json file with all the stats that can be used with the types graph visualization widget as describe in the output section.

## Output - hcls-output.rb

```
ruby hcls-output.rb


hcls-output.rb <hcld-dd_output.json> <prefixes.tsv>

# The script needs json output from the hcls-dd.rb script.
# Will produce TSV statistics table for each graphs and a graph.js file
# that can be used by the types graph viewer

```

We highly recommend to use prefixes for all the type and predicate to get a cleaner visualization output.

This will generate a directory rdf-stats/ with a subdirectory by graph (see ./example/rdf-stats/).

A graph directory contains 7 table files (TSV) :
* basic_metrics.tsv
* datatypes.tsv
* object-type_property.tsv
* properties.tsv
* subject-type_property.tsv
* type-type_property.tsv
* types.tsv

and 1 javascript file : type-type_property_graph.js

Replace the ./type-graphs-html/graph.js with the type-type_property_graph.js and open ./type-graphs-html/index.html file to visualize it.


## Example

See the ./example/ directory


## TODO

* Test and validate the RDF graph output
* Merge hcls-output.rb into hcls-dd.rb and add CLI options

