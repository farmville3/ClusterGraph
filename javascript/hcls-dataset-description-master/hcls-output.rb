#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# date:    	2016-06-01
# version: 	0.01

require 'json'
require 'digest'

usage = "
hcls-output.rb <hcld-dd_output.json> <prefixes.tsv>

# The script needs json output from the hcls-dd.rb script.
# Will produce TSV statistics table for each graphs and a graph.js file
# that can be used by the types graph viewer 

"

if ARGV.length < 1
  abort usage
end

prefixes = {}

if ARGV[1]
  File.open(ARGV[1]) do |f|
    while l = f.gets
      lA = l.chomp.split("\t")
      prefixes[lA[0]] = lA[1]
    end
  end
end

# ignore virtuoso graphs
virtuoso_graphs = [
  "http://www.openlinksw.com/schemas/virtrdf#",
  "http://www.w3.org/ns/ldp#",
  "http://URIQAREPLACEME/sparql",
  "http://URIQAREPLACEME/DAV/",
  "http://www.w3.org/2002/07/owl#",
  "b3sonto",
  "b3sifp",
  "urn:rules.skos",
  "http://www.openlinksw.com/schemas/oplweb#",
  "virtrdf-label",
  "facets",
  "_all"
]

hcls = JSON.parse(File.read(ARGV[0]))

# output directory
root = "rdf-stats/"
if ! Dir.exists? root
  Dir.mkdir(root)  
end


# process each graph
hcls.each do |k,v|
  next if virtuoso_graphs.include? "#{k}"

  outdir = k.gsub("http://","")
  puts outdir
  basic_metrics = {}

  if ! Dir.exists? (root + outdir)
    Dir.mkdir(root + outdir)
  end

  v.each do |stat,obj|

    case stat
    when /^c_/

      basic_metrics[stat] = obj

    when "types"

      outfile = File.open(root + outdir + "/types.tsv", "w")
      outfile.write("URI\tCount\tDistinct_count\n")
      outarray = []
      obj.each do |k2,v2|
        outarray << {
          uri: v2['uri'],
          distinct_count: v2['distinct_count'],
          count: v2['count'] }
      end
      outarray.sort_by{ |v2| -v2[:count]}.each do |v3|
        outfile.write("#{v3[:uri]}\t#{v3[:count]}\t#{v3[:distinct_count]}\n")
      end

    when "properties"

      outfile = File.open(root + outdir + "/properties.tsv", "w")
      outfile.write("URI\tCount\n")
      outarray = []
      obj.each do |k2,v2|
        outarray << {
          uri: v2['uri'],
          count: v2['count']
        }
      end
      outarray.sort_by{ |v2| -v2[:count]}.each do |v3|
        outfile.write("#{v3[:uri]}\t#{v3[:count]}\n")
      end

    when "datatypes"

      outfile = File.open(root + outdir + "/datatypes.tsv", "w")
      outfile.write("URI\tCount\tDistinct_count\n")
      outarray = []
      obj.each do |k2,v2|
        outarray << {
          uri: v2['uri'],
          count: v2['count'],
          distinct_count: v2['distinct_count']
        }
      end
      outarray.sort_by{ |v2| -v2[:count]}.each do |v3|
        outfile.write("#{v3[:uri]}\t#{v3[:count]}\t#{v3[:distinct_count]}\n")
      end

    when "subjecttype_property"

      outfile = File.open(root + outdir + "/subject-type_property.tsv", "w")
      outfile.write("Subject-Type\tProperty\tCount\tDistinct_count\n")
      outarray = []
      obj.each do |v2|
        outarray << {
          stype: v2['stype'],
          property: v2['property'],
          count: v2['count'],
          distinct_count: v2['distinct_count']
        }
      end
      outarray.sort_by{ |v2| -v2[:count]}.each do |v3|
        outfile.write("#{v3[:stype]}\t#{v3[:property]}\t#{v3[:count]}\t#{v3[:distinct_count]}\n")
      end

    when "objecttype_property"

      outfile = File.open(root + outdir + "/object-type_property.tsv", "w")
      outfile.write("Object-Type\tProperty\tCount\tDistinct_count\n")
      outarray = []
      obj.each do |v2|
        outarray << {
          otype: v2['otype'],
          property: v2['property'],
          count: v2['count'],
          distinct_count: v2['distinct_count']
        }
      end
      outarray.sort_by{ |v2| -v2[:count]}.each do |v3|
        outfile.write("#{v3[:otype]}\t#{v3[:property]}\t#{v3[:count]}\t#{v3[:distinct_count]}\n")
      end

    when "type_type_property"

      # node and edges for graph vis.js
      nodes = {}
      edges = []
      edges_it = 0

      outfile = File.open(root + outdir + "/type-type_property.tsv", "w")
      outfile.write("Subject-Type\tSubject-Type Count\tSubject-Type Distinct Count\tProperty\tObject-Type\tObject-Type Count\tObject-Type tDistinct_count\n")
      outarray = []

      obj.each do |v2|
        outarray << {
          stype: v2['stype'],
          property: v2['property'],
          s_count: v2['stype_count'],
          s_distinct_count: v2['stype_distinct_count'],
          otype: v2['otype'],
          o_count: v2['otype_count'],
          o_distinct_count: v2['otype_distinct_count']
        }

        # JSON Graph building for Vis.js
        weight = v2['stype_count']

        id = Digest::MD5.hexdigest(v2['stype'].to_s).to_i(16).to_s[0..12]
        if ! nodes.has_key? id or weight > nodes[id][:weight]
          label = v2['stype']
          prefixes.each do |prefix,subprefix|
            if label.include? prefix
              label.gsub!(prefix,subprefix)
              label.gsub!(":",":\\n")
            end
          end
          nodes[id] = {weight: weight, value: "{id: #{id}, weight: #{weight}, label: \'#{label}\'}"}
        end

        id2 = Digest::MD5.hexdigest(v2['otype'].to_s).to_i(16).to_s[0..12]
        if ! nodes.has_key? id2 or weight > nodes[id2][:weight]
          label = v2['otype']
          prefixes.each do |prefix,subprefix|
            if label.include? prefix
              label.gsub!(prefix,subprefix)
              label.gsub!(":",":\\n")
            end
          end
          nodes[id2] = {weight: weight, value: "{id: #{id2}, weight: #{weight}, label: \'#{label}\'}"}
        end

        label = v2['property']
        prefixes.each do |prefix,subprefix|
          if label.include? prefix
            label.gsub!(prefix, subprefix)
            label.gsub!(":",":\\n")
          end
        end
        edges << "{id: #{edges_it}, from: #{id}, to: #{id2}, weight: #{weight}, label: \'#{label}\'}"
        edges_it += 1
      end

      # stats file
      outarray.sort_by{ |v2| -v2[:s_count]}.each do |v3|
        outfile.write("#{v3[:stype]}\t#{v3[:s_count]}\t#{v3[:s_distinct_count]}\t#{v3[:property]}\t#{v3[:otype]}\t#{v3[:o_count]}\t#{v3[:o_distinct_count]}\n")
      end
      outfile.close

      outjsongraph = File.open(root + outdir + "/type-type_property_graph.js", "w")
      # JSON graph
      outjsongraph.write("var nodes_array = [\n")
      nodes.each do |k, node|
        outjsongraph.write(node[:value]+",\n")
      end
      outjsongraph.write("]\n\n")
      outjsongraph.write("var edges_array = [\n")
      edges.each do |edge|
        outjsongraph.write(edge+",\n")
      end
      outjsongraph.write("]")
      outjsongraph.close
    end

  end

  outfile = File.open(root + outdir + "/basic_metrics.tsv", "w")
  basic_metrics.sort_by {|_key, _value| -_value}.each do |m,n|
    outfile.write("#{m}\t#{n}\n")
  end

end
