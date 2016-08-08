#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# date:    	2016-03-07
# version: 	0.01

# HCLS Dataset Description queryer (statistics queries)
# https://www.w3.org/TR/hcls-dataset/#s6_6

# maximize the decoupling of queries (by graph/type) to improve their viability & performance

require 'rdf'
require 'rdf/turtle'
require 'net/http'
require 'json'
require 'digest'

class EndpointStatistics

  def initialize cli_arg

    # create @arg, @uri, @sparql_url, @http
    parse_cmd_arg(cli_arg)

    @endpoint_stats = {}
    @void_graph = RDF::Graph.new(format: :ttl)

    list_graph

    logs = {}
    if @arg[:basic_metrics]
      logs = compute_basic_metrics
    else
      logs = compute_all_metrics
    end

    puts "#Logs.."
    logs.each do |k,v|
      puts "#{k} - #{v}"
    end

  end

  # parse command line arguments
  def parse_cmd_arg arguments

    usage = "
    hcls-dd.rb [options] <SPARQL url>

    Input:
       --uri               <endpoint uri>
       --sparql            <target sparql uri>
       --default-graph     <default graph uri>
       --label             <,\'endpoint label\'> [default: my_endpoint]
       --description       <\'endpoint description\'> [default: My Endpoint Statistics]

    Output:
       --output_json       <outfile.json> save hash as json file
       --output_rdf        <outfile.ttl>

    Options:
       --basic_metrics     to compute only the basic metrics

       --help              print this
    "

    if arguments.length < 1
      abort usage
    end

    # Set the HTTP agent
    @sparql_url = arguments[-1]
    abort usage if @sparql_url[0..1] == "--"
    abort "Need http:// to SPARQL URL !" if ! @sparql_url.include? "http://"
    url_split = @sparql_url.split("/")
    @http = ""
    if url_split[2].include?(":")
      @http = Net::HTTP.new(url_split[2].split(":")[0],url_split[2].split(":")[1])
    else
      @http = Net::HTTP.new(url_split[2])
    end

    @http.read_timeout = 360000000 # seconds
    @http.open_timeout = 360000000 # seconds

    @uri = "http://" + @sparql_url.gsub("http://","").split(":")[0]
    @sparql_uri = @uri+"/sparql"
    @label = "my_endpoint"

    @arg = {}
    @arg[:basic_metrics] = false

    # reading opts
    for i in 0..arguments.length-2
      next if ARGV[i][0] != "-"
      key = ARGV[i].gsub("-","").downcase
      if key == "uri"
        @uri = ARGV[i+1]
      elsif key == "sparql"
        @sparql_uri = ARGV[i+1]
      elsif key == "label"
        @label = ARGV[i+1]
      elsif key == "basic_metrics"
        @arg[:basic_metrics] = true
      # elsif key == "output_html"
      #   @arg[:html] = ARGV[i+1]
      elsif key == "help"
        abort usage
      else
        @arg[key.to_sym] = ARGV[i+1]
      end
    end

  end

  # list graph all graphs
  def list_graph
    puts "#List Graphs"
    graphs = []
    q = "SELECT DISTINCT $g WHERE { GRAPH $g { $s $p $o . }}"
    new_uri = @sparql_url+"?query=#{q}&format=json"
    req = @http.request_get(URI(new_uri))
    res = JSON.parse(req.body)
    res['results']['bindings'].each do |tpl|
      graphs << tpl['g']['value']
    end
    graphs.each do |g|
      @endpoint_stats[g] = {}
      puts "Graph\t#{g}"
    end
    @endpoint_stats['_all'] = {}
    graphs
  end

  # compute all basic metrics
  def compute_basic_metrics

    logs = {}

    all_metrics = [
      "sum_triples",
      "sum_subjects",
      "sum_properties",
      "sum_objects",
      "sum_entities",
      "sum_literals",
      "sum_types" ]

    try_count = 5

    all_metrics.each do |metric|
      i = 0

      begin

        puts "# Computing #{metric}.."
        send(metric)
        logs[metric] = "success\t#{i+1}"

      rescue => e

        i+=1
        retry if i<try_count
        puts "FAILED.. aborting this metric ; #{metric}\n#{e}"
        logs[metric] = "#{e}"

      end

    end

    logs
    
  end


  # compute all metrics
  def compute_all_metrics

    logs = {}

    all_metrics = [
      "sum_triples",
      "sum_subjects",
      "sum_properties",
      "sum_objects",
      "sum_entities",
      "sum_literals",
      "sum_types",
      "type_count",
      "property_count",
      "datatype_count",
      "subjecttype_property_count",
      "objecttype_property_count",
      "type_type_count",
      "dataset_dataset_count"]

    try_count = 5

    all_metrics.each do |metric|
      i = 0

      begin

        puts "# Computing #{metric}.."
        send(metric)
        logs[metric] = "success\t#{i+1}"

      rescue => e

        i+=1
        retry if i<try_count
        puts "FAILED.. aborting this metric ; #{metric}\n#{e}"
        logs[metric] = "#{e}"

      end

    end

    logs
    
  end

  # sum total number of triples [by graph]
  def sum_triples
    puts "#Number of Triples"
    nb_of_triples = 0
    @endpoint_stats.each_key do |k|
      next if k == "_all"
      q = "SELECT COUNT(*) as $nb WHERE { GRAPH <#{URI.encode(k)}> { $s $p $o . }}"
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        @endpoint_stats[k]['c_triples'] = tpl['nb']['value'].to_i
        nb_of_triples += tpl['nb']['value'].to_i
        puts "Triples\t[#{k}]\t#{tpl['nb']['value']}"
      end
    end  
    @endpoint_stats['_all']['c_triples'] = nb_of_triples
    nb_of_triples
  end

  # sum number of distinct literals [by graph]
  def sum_literals
    puts "#Number of Literals"
    nb_of_literals = 0
    @endpoint_stats.each_key do |k|
      next if k == "_all"
      g = "<#{URI.encode(k)}>"
      q = "SELECT COUNT (DISTINCT $o) as $nb WHERE { GRAPH #{g} { $s $p $o . FILTER (isLiteral($o))}}"
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        @endpoint_stats[k]['c_literals'] = tpl['nb']['value'].to_i
        nb_of_literals += tpl['nb']['value'].to_i
        puts "Literals\t[#{k}]\t#{tpl['nb']['value']}"
      end
    end
    @endpoint_stats['_all']['c_literals'] = nb_of_literals
    nb_of_literals
  end

  # sum number of distinct objects [by graph]
  def sum_objects
    puts "#Number of Objects"
    nb_of_objects = 0
    @endpoint_stats.each_key do |k|
      next if k == "_all"
      g = "<#{URI.encode(k)}>"
      q = "SELECT COUNT (DISTINCT $o) as $nb WHERE { GRAPH #{g} { $s $p $o . FILTER (!isLiteral($o))}}"
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        @endpoint_stats[k]['c_objects'] = tpl['nb']['value'].to_i
        nb_of_objects += tpl['nb']['value'].to_i
        puts "Objects\t[#{k}]\t#{tpl['nb']['value']}"
      end
    end
    @endpoint_stats['_all']['c_objects'] = nb_of_objects
    nb_of_objects
  end

  # sum distinct entities
  def sum_entities
    puts "#Number of Entities"
    last_count = 0
    @endpoint_stats.each_key do |k|
      next if k == "_all"
      g = "<#{URI.encode(k)}>"
      q = "SELECT COUNT (DISTINCT $s) as $nb WHERE { GRAPH #{g} { $s a [] . }}"
      nb_entities = -1
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        nb_entities = tpl['nb']['value'].to_i
      end
      @endpoint_stats[k]['c_entities'] = nb_entities
      puts "Entities #{k} : " + "#{@endpoint_stats[k]['c_entities']}"
    end
  end

  # sum distinct subjects
  def sum_subjects
    puts "#Number of Subjects"
    last_count = 0
    @endpoint_stats.each_key do |k|
      next if k == "_all"
      g = "<#{URI.encode(k)}>"
      q = "SELECT COUNT (DISTINCT $s) as $nb WHERE { GRAPH #{g} { $s $p $o . }}"
      nb_subjects = -1
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        nb_subjects = tpl['nb']['value'].to_i
      end
      @endpoint_stats[k]['c_subjects'] = nb_subjects
      puts "Subjects #{k} : " + "#{@endpoint_stats[k]['c_subjects']}"
    end
  end

  # sum distinct property
  def sum_properties
    puts "#Number of Properties"
    last_count = 0
    @endpoint_stats.each_key do |k|
      next if k == "_all"
      g = "<#{URI.encode(k)}>"
      q = "SELECT COUNT (DISTINCT $p) as $nb WHERE { GRAPH #{g} { $s $p $o . }}"
      nb_properties = -1
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        nb_properties = tpl['nb']['value'].to_i
      end
      @endpoint_stats[k]['c_properties'] = nb_properties
      puts "Properties #{k} : " + "#{@endpoint_stats[k]['c_properties']}"
    end
  end

  # sum number of distinct types / class [by graph]
  def sum_types
    puts "#Number of Types"
    @endpoint_stats.each_key do |k|
      @endpoint_stats[k]['types'] = {}
      next if k == "_all"
      g = "<#{URI.encode(k)}>"
      q = "SELECT DISTINCT $o WHERE { GRAPH #{g} { $s a $o .}}"
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        @endpoint_stats[k]['types'][tpl['o']['value'].to_s] = {
          'uri' => tpl['o']['value'],
          'count' => 0,
          'distinct_count' => 0,
          'label' => "" }
      end
      puts "Types\t[#{k}]\t#{@endpoint_stats[k]['types'].length}"
      @endpoint_stats[k]['c_types'] = @endpoint_stats[k]['types'].length
    end
    types = @endpoint_stats.collect {|k,v| v['types'].keys }.flatten(1).uniq
    @endpoint_stats['_all']['types'] = types
    @endpoint_stats['_all']['c_types'] = types.length
    types.length
  end

  # count the number of entities by type
  def type_count
    puts "#Types count"
    @endpoint_stats.each_key do |k|
      next if k=="_all"
      g = "<#{URI.encode(k)}>"
      @endpoint_stats[k]['types'].each_key do |t|
        type = "<#{URI.encode(t)}>"
        q = "SELECT distinct $type $n $dn (str($label) AS $slabel) { GRAPH #{g} { SELECT $type (COUNT(DISTINCT $s) AS $dn) (COUNT($s) AS $n) {$s a $type . FILTER ($type=#{type})}} OPTIONAL {$type rdfs:label $label}}"
        new_uri = @sparql_url+"?query=#{q}&format=json"
        req = @http.request_get(URI(new_uri))
        res = JSON.parse(req.body)
        res["results"]["bindings"].each do |tpl|
          typek = tpl['type']['value']
          @endpoint_stats[k]['types'][typek]['count'] = tpl['n']['value'].to_i
          @endpoint_stats[k]['types'][typek]['distinct_count'] = tpl['dn']['value'].to_i
          if tpl.has_key? 'slabel'
            @endpoint_stats[k]['types'][typek]['label'] = tpl['slabel']['value']
          else
            @endpoint_stats[k]['types'][typek]['label'] = ""
          end
        end
      end
      puts "\t[#{k}]"
    end
  end

  # property count
  def property_count
    puts "#Properties count"
    @endpoint_stats.each_key do |k|
      next if k=="_all"
      @endpoint_stats[k]['properties'] = {}
      g = "<#{URI.encode(k)}>"
      q = "SELECT DISTINCT $p (str($plabel) AS $plabel) $n { GRAPH #{g} { SELECT $p (COUNT($p) AS $n) {$s $p $o} GROUP BY $p } OPTIONAL {$p rdfs:label $plabel}}"
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        prop = tpl['p']['value']
        @endpoint_stats[k]['properties'][prop] = {
          'uri' => prop,
          'count' => tpl['n']['value'].to_i,
          'label' => ""
        }
        if tpl.has_key? 'plabel'
          @endpoint_stats[k]['properties'][prop]['label'] = tpl['plabel']['value']
        end
      end
      puts "\t[#{k}]"
    end
  end

  # property count (object property count)
  def object_property_count
    puts "#Object Properties count"
    @endpoint_stats.each_key do |k|
      next if k=="_all"
      @endpoint_stats[k]['object_properties'] = {}
      g = "<#{URI.encode(k)}>"
      q = "SELECT DISTINCT $p (str($label) AS $plabel) ($n AS $n) ($dn AS $dn){ GRAPH #{g} { SELECT $p (COUNT($o) AS $n) (COUNT(DISTINCT $o) AS $dn) { $s $p $o FILTER (!isLiteral($o)) } GROUP BY $p} OPTIONAL {$p rdfs:label $label}}"
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        prop = tpl['p']['value']
        @endpoint_stats[k]['object_properties'][prop] = {
          'uri' => prop,
          'count' => tpl['n']['value'].to_i,
          'distinct_count' => tpl['dn']['value'].to_i,
          'label' => ""
        }
        if tpl.has_key? 'plabel'
          @endpoint_stats[k]['object_properties'][prop]['label'] = tpl['plabel']['value']
        end
      end
      puts "\t[#{k}]"
    end
  end


  # datatype count
  def datatype_count
    puts "#Datatypes count"
    @endpoint_stats.each_key do |k|
      next if k=="_all"
      @endpoint_stats[k]['datatypes'] = {}
      g = "<#{URI.encode(k)}>"
      q = "SELECT DISTINCT $p (str($label) AS $plabel) ($n AS $n) ($dn AS $dn){ GRAPH #{g} { SELECT $p (COUNT($o) AS $n) (COUNT(DISTINCT $o) AS $dn){$s $p $o . FILTER(isLiteral($o))} GROUP BY $p} OPTIONAL {$p rdfs:label $label}}"
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        prop = tpl['p']['value']
        @endpoint_stats[k]['datatypes'][prop] = {
          'uri' => prop,
          'count' => tpl['n']['value'].to_i,
          'distinct_count' => tpl['dn']['value'].to_i,
          'label' => ""
        }
        if tpl.has_key? 'plabel'
          @endpoint_stats[k]['datatypes'][prop]['label'] = tpl['plabel']['value']
        end
      end
      puts "\t[#{k}]"
    end
  end

  # subjecttype + property count
  def subjecttype_property_count
    puts "#Subject-Type_Property count"
    @endpoint_stats.each_key do |k|
      next if k=="_all"
      @endpoint_stats[k]['subjecttype_property'] = []
      g = "<#{URI.encode(k)}>"
      q = "SELECT DISTINCT $p (str($plabel) AS $plabel) $stype (str($stype_label) AS $stype_label) ($n AS $n) ($dn AS $dn){ GRAPH #{g} { SELECT $p $stype (COUNT($s) AS $n) (COUNT(DISTINCT $s) AS $dn) { $s $p $o . $s a $stype . } GROUP BY $p $stype } OPTIONAL {$p rdfs:label $plabel} OPTIONAL {$stype rdfs:label $stype_label}}"
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        stype_prop = {
          'stype' => tpl['stype']['value'],
          'property' => tpl['p']['value'],
          'count' => tpl['n']['value'].to_i,
          'distinct_count' => tpl['dn']['value'].to_i,
          'stype_label' => "",
          'p_label' => ""
        }
        stype_prop['stype_label'] = tpl['stype_label']['value'] if tpl.has_key? 'stype_label'
        stype_prop['p_label'] = tpl['p_label']['value'] if tpl.has_key? 'p_label'
        @endpoint_stats[k]['subjecttype_property'] << stype_prop
      end
      puts "\t[#{k}]"
    end
  end

  # objecttype + property count
  def objecttype_property_count
    puts "#Object-Type_Property count"
    @endpoint_stats.each_key do |k|
      next if k=="_all"
      @endpoint_stats[k]['objecttype_property'] = []
      g = "<#{URI.encode(k)}>"
      q = "SELECT distinct $p (str($plabel) AS $plabel) $otype (str($otype_label) AS $otype_label) ($n AS $n) ($dn AS $dn) { GRAPH #{g} { SELECT $p $otype (COUNT($o) AS $n) (COUNT(DISTINCT $o) AS $dn) { $s $p $o . $o a $otype . } GROUP BY $p $otype } OPTIONAL {$p rdfs:label $plabel} OPTIONAL {$otype rdfs:label $otype_label}}"
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        otype_prop = {
          'otype' => tpl['otype']['value'],
          'property' => tpl['p']['value'],
          'count' => tpl['n']['value'].to_i,
          'distinct_count' => tpl['dn']['value'].to_i,
          'otype_label' => "",
          'p_label' => ""
        }      
        otype_prop['otype_label'] = tpl['otype_label']['value'] if tpl.has_key? 'otype_label'
        otype_prop['p_label'] = tpl['p_label']['value'] if tpl.has_key? 'p_label'
        @endpoint_stats[k]['objecttype_property'] << otype_prop
      end
      puts "\t[#{k}]"
    end
  end


  # type - type count
  def type_type_count
    puts "#Type-Type Property count"
    @endpoint_stats.each_key do |k|
      next if k=="_all"
      @endpoint_stats[k]['type_type_property'] = []
      predicates = []
      g = "<#{URI.encode(k)}>"
      q1 = "SELECT $p { GRAPH #{g} { $s $p $o FILTER (!isLiteral($o)) }} GROUP BY $p"
      new_uri = @sparql_url+"?query=#{q1}&format=json"
      req1 = @http.request_get(URI(new_uri))
      res1 = JSON.parse(req1.body)
      res1["results"]["bindings"].each do |tpl|
        predicates << tpl['p']['value']
      end
      predicates.each do |p|
        type_type_pred = {}
        pred = "<#{URI.encode(p)}>"
        q2 = "SELECT DISTINCT $stype (str($stype_label) AS $stype_label) ($sn AS $sn) ($dsn AS $dsn) $p (str($plabel) AS $plabel) $otype (str($otype_label) AS $otype_label) ($on AS $on) ($don AS $don) { GRAPH #{g} { SELECT distinct $stype $p $otype (COUNT($s) AS $sn) (COUNT(DISTINCT $s) AS $dsn) (COUNT($o) AS $on) (COUNT(DISTINCT $o) AS $don){ $s $p $o . $s a $stype . $o a $otype . FILTER($p = #{pred})} GROUP BY $p $stype $otype } OPTIONAL {$stype rdfs:label $stype_label} OPTIONAL {$p rdfs:label $plabel} OPTIONAL {$otype rdfs:label $otype_label}}"
        new_uri = @sparql_url+"?query=#{q2}&format=json"
        req2 = @http.request_get(URI(new_uri))
        res2 = JSON.parse(req2.body)
        res2["results"]["bindings"].each do |tpl|
          type_type  = {
            'stype' => tpl['stype']['value'],
            'stype_count'=> tpl['sn']['value'].to_i,
            'stype_distinct_count' => tpl['dsn']['value'].to_i,
            'stype_label' => "",
            'property' => p,
            'otype' => tpl['otype']['value'],
            'otype_count'=> tpl['on']['value'].to_i,
            'otype_distinct_count' => tpl['don']['value'].to_i,
            'otype_label' => ""
          }
          type_type['stype_label'] = tpl['stype_label']['value'] if tpl.has_key? ['stype_label']
          type_type['otype_label'] = tpl['otype_label']['value'] if tpl.has_key? ['otype_label']      
          @endpoint_stats[k]['type_type_property'] << type_type
        end
      end
      puts "\t[#{k}]"
    end
  end

  # dataset - dataset
  def dataset_dataset_count
    puts "#Dataset-Dataset Property count"
    @endpoint_stats.each_key do |k|
      next if k=="_all"
      @endpoint_stats[k]['dataset_dataset_property'] = []
      g = "<#{URI.encode(k)}>"
      dataset_dataset = {}
      q = "SELECT DISTINCT $p $stype $otype (COUNT($s) AS $n) { GRAPH #{g} { $s $p $o . $s a $stype . $o a $otype . FILTER regex ($stype, 'vocabulary:Resource') FILTER regex ($otype, 'vocabulary:Resource') FILTER ($stype != $otype)}}"
      new_uri = @sparql_url+"?query=#{q}&format=json"
      req = @http.request_get(URI(new_uri))
      res = JSON.parse(req.body)
      res["results"]["bindings"].each do |tpl|
        dataset_dataset = {
          'stype' => tpl['stype']['value'],
          'otype' => tpl['otype']['value'],
          'property' => tpl['p']['value'],
          'count' => tpl['n']['value'].to_i
        }
        @endpoint_stats[k]['dataset_dataset_property'] << dataset_dataset
      end
      puts "\t[#{k}]"
    end
  end


  # create a classpartition for VOID/HCLS
  def class_partition_bgp graph_void, void_p, type, void_c, value
    triples = []
    res = Digest::SHA256.hexdigest "#{type}+#{value}"
    uri = RDF::URI("#{@uri}/dataset_resource:#{res[-16..-1]}")
    triples << RDF::Statement.new(
      graph_void,
      RDF::URI("http://rdfs.org/ns/void#classPartition"),
      uri
    )
    triples << RDF::Statement.new(
      uri,
      RDF::URI("http://rdfs.org/ns/void##{void_p}"),
      RDF::URI(type)
    )
    triples << RDF::Statement.new(
      uri,
      RDF::URI("http://rdfs.org/ns/void##{void_c}"),
      RDF::Literal(value)
    )
    triples
  end


  # create rdf graph of the statistics (void+hcls_stats)
  # See -> Using VoID with the SPARQL Service Description Vocabulary (https://www.w3.org/TR/void/#sparql-sd)
  # See -> Dataset Descriptions: HCLS Community Profile  (https://www.w3.org/TR/hcls-dataset/#s6_6)
  def create_statistics_graph

    puts "#Ouput statistics graph"
    sparql_sd = "http://www.w3.org/ns/sparql-service-description#"
    void = "http://rdfs.org/ns/void#"
    rdfs = "http://www.w3.org/2000/01/rdf-schema#"
    void_ext  = "http://ldf.fi/void-ext#"

    dataset_bnode = RDF::URI(@uri+"/DatasetDescription")

    sparql_service = []

    sparql_service << RDF::Statement.new(RDF::URI(@uri),RDF.type,RDF::URI(sparql_sd+"Service"))
    sparql_service << RDF::Statement.new(RDF::URI(@uri),RDF::URI(sparql_sd+"url"),RDF::URI(@sparql_uri))
    sparql_service << RDF::Statement.new(RDF::URI(@uri),RDF::URI(sparql_sd+"defaultDatasetDescription"),dataset_bnode)
    sparql_service << RDF::Statement.new(dataset_bnode,RDF.type,RDF::URI(sparql_sd+"Dataset"))

    sparql_service.each do |t|
      if t.valid?
        @void_graph << t
      end
    end

    # Each graph -> create a named graph
    @endpoint_stats.each do |k,v|
      triples = []

      puts "#{k}"

      graph_void_uri = "#{k}/VOID"
      graph_void = RDF::URI("#{k}/VOID")
      graph = RDF::URI("#{k}")

      triples << RDF::Statement.new(dataset_bnode,RDF::URI(sparql_sd+"namedGraph"),graph)
      triples << RDF::Statement.new(graph,RDF::URI(sparql_sd+"name"), graph)
      triples << RDF::Statement.new(graph,RDF::URI(sparql_sd+"graph"),graph_void)
      triples << RDF::Statement.new(graph_void,RDF.type,RDF::URI(sparql_sd+"Graph"))
      triples << RDF::Statement.new(graph_void,RDF.type,RDF::URI(void+"Dataset"))

      # Basic metric count
      # 6.6.1.1 To specify the number of triples in the dataset
      triples << RDF::Statement.new(graph_void,RDF::URI(void+"triples"),RDF::Literal(v['c_triples']))
      # 6.6.1.2 To specify the number of unique, typed entities in the dataset
      triples << RDF::Statement.new(graph_void,RDF::URI(void+"entities"),RDF::Literal(v['c_entities']))
      # 6.6.1.3 To specify the number of unique subjects in the dataset
      triples << RDF::Statement.new(graph_void,RDF::URI(void+"distinctSubjects"),RDF::Literal(v['c_subjects']))
      # 6.6.1.4 To specify the number of unique properties in the dataset
      triples << RDF::Statement.new(graph_void,RDF::URI(void+"properties"),RDF::Literal(v['c_properties']))
      # 6.6.1.5 To specify the number of unique objects in the dataset
      triples << RDF::Statement.new(graph_void,RDF::URI(void+"distinctObjects"),RDF::Literal(v['c_objects']))

      # Number of Classes / Types
      # 6.6.1.6 To specify the number of unique classes in the dataset
      res = Digest::SHA256.hexdigest "#{k}+classes_count+#{v['c_types']}"
      res_uri = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res[-16..-1]}")
      triples << RDF::Statement.new(graph_void,RDF::URI(void+"classPartition"),res_uri)
      triples << RDF::Statement.new(res_uri,RDF::URI(void+"class"),RDF::URI(rdfs+"Class"))
      triples << RDF::Statement.new(res_uri,RDF::URI(void+"distinctSubjects"),RDF::Literal(v['c_types']))

      # Number of Literals
      # 6.6.1.7 To specify the number of unique literals in the dataset
      res = Digest::SHA256.hexdigest "#{k}+literals_count+#{v['c_literals']}"
      res_uri = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res[-16..-1]}")
      triples << RDF::Statement.new(graph_void,RDF::URI(void+"classPartition"),res_uri)
      triples << RDF::Statement.new(res_uri,RDF::URI(void+"class"),RDF::URI(rdfs+"Literals"))
      triples << RDF::Statement.new(res_uri,RDF::URI(void+"distinctSubjects"),RDF::Literal(v['c_literals']))

      # Number of Graph (TODO)
      # 6.6.1.8 To specify the number of graphs in the dataset

      # Count by obj type
      # 6.6.2.1 To specify the classes and the number of their instances in the dataset
      if v.has_key? 'types'
        v['types'].each do |type_k, type_v|
          next if type_v.nil?
          res = Digest::SHA256.hexdigest "#{k}+#{type_k}+classes_count+#{type_v['distinct_count']}"
          res_uri = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res[-16..-1]}")
          triples << RDF::Statement.new(graph_void,RDF::URI(void+"classPartition"),res_uri)
          triples << RDF::Statement.new(res_uri,RDF::URI(void+"class"),RDF::URI(type_v['uri']))
          triples << RDF::Statement.new(res_uri,RDF::URI(void+"distinctSubjects"),RDF::Literal(type_v['distinct_count']))
          triples << RDF::Statement.new(res_uri,RDF::URI(void+"triples"),RDF::Literal(type_v['count']))
        end
      end

      # Count by property
      # 6.6.2.2 To specify the properties and their occurrence in the dataset
      if v.has_key? 'properties'
        v['properties'].each do |prop_k, prop_v|
          next if prop_v.nil?
          res = Digest::SHA256.hexdigest "#{k}+#{prop_k}+triples_count+#{prop_v['count']}"
          res_uri = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res[-16..-1]}")
          triples << RDF::Statement.new(graph_void,RDF::URI(void+"propertyPartition"),res_uri)
          triples << RDF::Statement.new(res_uri,RDF::URI(void+"property"),RDF::URI(prop_v['uri']))
          triples << RDF::Statement.new(res_uri,RDF::URI(void+"triples"),RDF::Literal(prop_v['count']))
        end
      end

      # Count by datatypes
      # TODO which one is this
      # if v.has_key? 'datatypes'
      #   v['datatypes'].each do |datatype_k, datatype_v|
      #     next if datatype_v.nil?
      #     res = Digest::SHA256.hexdigest "#{k}+#{datatype_k}+triples_count+#{v['c_datatypes']}"
      #     res_uri = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res[-16..-1]}")
      #     triples << RDF::Statement.new(graph_void,RDF::URI(void+"classPartition"),res_uri)
      #     triples << RDF::Statement.new(res_uri,RDF::URI(void+"property"),RDF::URI(datatype_v['uri']))
      #     triples << RDF::Statement.new(res_uri,RDF::URI(void+"triples"),RDF::Literal(datatype_v['count']))
      #     triples << RDF::Statement.new(res_uri,RDF::URI(void+"distinctSubjects"),RDF::Literal(datatype_v['distinct_count']))
      #   end
      # end


      # Count subject type associated to properties
      # 6.6.2.3 To specify the property, the number of unique typed subjects, and number of triples linked to a property in the dataset
      if v.has_key? 'subjecttype_property'
        # stype_prop = {
        #   'stype' => tpl['stype']['value'],
        #   'property' => tpl['p']['value'],
        #   'count' => tpl['n']['value'].to_i,
        #   'distinct_count' => tpl['dn']['value'].to_i,
        #   'stype_label' => "",
        #   'p_label' => ""
        # }
        v['subjecttype_property'].each do |stype_prop_v|
          next if stype_prop_v.nil?
          res = Digest::SHA256.hexdigest "#{k}+#{stype_prop_v['property']}+#{stype_prop_v['stype']}+triples_count"
          res_uri = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res[-16..-1]}")
          triples << RDF::Statement.new(graph_void,RDF::URI(void+"propertyPartition"),res_uri)
          triples << RDF::Statement.new(res_uri,RDF::URI(void+"property"),RDF::URI(stype_prop_v['property']))
          triples << RDF::Statement.new(res_uri,RDF::URI(void+"triples"),RDF::Literal(stype_prop_v['count']))
          res2 = Digest::SHA256.hexdigest "#{k}+#{stype_prop_v['property']}+#{stype_prop_v['stype']}+distinctSubject"
          res_uri2 = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res2[-16..-1]}")
          triples << RDF::Statement.new(res_uri, RDF::URI(void+"classPartition"),res_uri2)
          triples << RDF::Statement.new(res_uri2,RDF::URI(void+"class"),RDF::URI(stype_prop_v['stype']))
          triples <<  RDF::Statement.new(res_uri2,RDF::URI(void+"distinctSubjects"),RDF::Literal(stype_prop_v['distinct_count']))
        end
      end

      # # Count object types associated to properties
      # # 6.6.2.4 To specify the number of unique typed objects linked to a property in the dataset:
      if v.has_key? 'objecttype_property'
        # otype_prop = {
        #   'otype' => tpl['otype']['value'],
        #   'property' => tpl['p']['value'],
        #   'count' => tpl['n']['value'].to_i,
        #   'distinct_count' => tpl['dn']['value'].to_i,
        #   'otype_label' => "",
        #   'p_label' => ""
        # }
        v['objecttype_property'].each do |otype_prop_v|
          next if otype_prop_v.nil?
          res = Digest::SHA256.hexdigest "#{k}+#{otype_prop_v['property']}+#{otype_prop_v['otype']}+triples_count"
          res_uri = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res[-16..-1]}")
          triples << RDF::Statement.new(graph_void,RDF::URI(void+"propertyPartition"),res_uri)
          triples << RDF::Statement.new(res_uri,RDF::URI(void+"property"),RDF::URI(otype_prop_v['property']))
          triples << RDF::Statement.new(res_uri,RDF::URI(void+"triples"),RDF::Literal(otype_prop_v['count']))
          res2 = Digest::SHA256.hexdigest "#{k}+#{otype_prop_v['property']}+#{otype_prop_v['otype']}+distinctSubject"
          res_uri2 = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res2[-16..-1]}")
          triples << RDF::Statement.new(res_uri, RDF::URI(void_ext+"objectClassPartition"),res_uri2)
          triples << RDF::Statement.new(res_uri2,RDF::URI(void+"class"),RDF::URI(otype_prop_v['otype']))
          triples <<  RDF::Statement.new(res_uri2,RDF::URI(void+"distinctObjects"),RDF::Literal(otype_prop_v['distinct_count']))
        end
      end

      # Count Literals associated to property TODO
      # 6.6.2.5 To specify the triples and number of unique literals that are related to a property in the dataset

      # Count unique subject + unique subject associated to properties
      # 6.6.2.6 To specify the number of unique subject types that are linked to unique object types in the dataset
      if v.has_key? 'type_type'
        v['type_type'].each do |type_type_v|
          # type_type  = {
          #   'stype' => tpl['stype']['value'],
          #   'stype_count'=> tpl['sn']['value'].to_i,
          #   'stype_distinct_count' => tpl['dsn']['value'].to_i,
          #   'stype_label' => "",
          #   'property' => p,
          #   'otype' => tpl['otype']['value'],
          #   'otype_count'=> tpl['on']['value'].to_i,
          #   'otype_distinct_count' => tpl['don']['value'].to_i,
          #   'otype_label' => ""
          # }
          next if type_type_v.nil?
          res = Digest::SHA256.hexdigest "#{k}+#{type_type_v['property']}+#{type_type_v['otype']}+triples_count"
          res_uri = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res[-16..-1]}")
          triples << RDF::Statement.new(graph_void,RDF::URI(void+"propertyPartition"),res_uri)
          triples << RDF::Statement.new(res_uri,RDF::URI(void+"property"),RDF::URI(type_type_v['property']))
          triples << RDF::Statement.new(res_uri,RDF::URI(void+"triples"),RDF::Literal(type_type_v['count']))
          res2 = Digest::SHA256.hexdigest "#{k}+#{type_type_v['property']}+#{type_type_v['otype']}+distinctSubject"
          res_uri2 = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res2[-16..-1]}")
          triples << RDF::Statement.new(res_uri, RDF::URI(void_ext+"objectClassPartition"),res_uri2)
          triples << RDF::Statement.new(res_uri2,RDF::URI(void+"class"),RDF::URI(type_type_v['otype']))
          triples <<  RDF::Statement.new(res_uri2,RDF::URI(void+"distinctObjects"),RDF::Literal(type_type_v['distinct_count']))
          res_uri3 = RDF::URI("#{@uri}/#{@label}.dataset_resource:#{res2[-16..-1]}")
          triples << RDF::Statement.new(res_uri, RDF::URI(void_ext+"objectClassPartition"),res_uri3)
          triples << RDF::Statement.new(res_uri3,RDF::URI(void+"class"),RDF::URI(type_type_v['otype']))
          triples <<  RDF::Statement.new(res_uri3,RDF::URI(void+"distinctObjects"),RDF::Literal(type_type_v['distinct_count']))
          
        end
      end

      triples.each do |t|
        if t.valid?
          @void_graph << t
        end
      end

    end

    # Write graph to disk
    filename = "output.ttl"
    if @arg.has_key? :output_rdf
      filename = @arg[:output_rdf]
    end

    fout = File.open(filename,"w")
    fout.close
    
    RDF::Turtle::Writer.open(filename,
                             prefixes: {
                               void: void,
                               void_ext: void_ext,
                               sparql_sd: sparql_sd,
                               rdfs: rdfs,
                               @label.to_sym => "#{@uri}/",
                             }) do |writer|

      writer << @void_graph
    end

  end

  # save json hash file
  def create_json_hash
    #save json
    output_file = "output.json"
    if @arg.has_key? :output_json
      output_file = @arg[:output_json]
    end
    File.open("#{output_file}","w") do |fout|
      fout.write(JSON.pretty_generate(@endpoint_stats))
    end
  end

end


# # Main # #

endpoint_stats = EndpointStatistics.new(ARGV)

# logs = endpoint_stats.compute_all_metrics

endpoint_stats.create_statistics_graph
endpoint_stats.create_json_hash
