var nodes, edges, network, neighbor_nodes, selected_node;

function addNode(node) {
  try {
    nodes.add(node);
  }
  catch (err) {
    alert(err);
  }
}

function updateNode(node) {
  try {
    nodes.update(node);
  }
  catch (err) {
    alert(err);
  }
}

function removeNode(id) {
  try {
    nodes.remove({id: id});
    index = _.findIndex(nodes_array,{id: id})
    nodes_array.splice(index)
  }
  catch (err) {
    alert(err);
  }
}

function addEdge(edge) {
  try {
    edges.add(edge);
  }
  catch (err) {
    alert(err);
  }
}

function updateEdge(edge) {
  try {
    edges.update(edge);
  }
  catch (err) {
    alert(err);
  }
}

function removeEdge(id) {
  try {
    edges.remove({id: id});
    index = _.findIndex(edges_array,{id: id})
    edges_array.splice(index)
  }
  catch (err) {
    alert(err);
  }
}

// Group functions
function removeNodes(group) {
  for (var i=0; i < nodes_array.length; i++) {
    if (nodes_array[i].group === group) {
      removeNode(nodes_array[i].id)
    }
  }
  // console.log("deleting nodes from group " + group)
}
function addNodes(group) {
  for (var i=0; i < nodes_array.length; i++) {
    if (nodes_array[i].group === group) {
      addNode(nodes_array[i])
    }
  }
  // console.log("adding nodes from group " + group)
}

// Selection functions
function removeSelectedNode(id) {
  removeNode(id);
  $('#selection').html("");
}

function nodeInfo(params) {

  var node_label = "";
  var incoming_edges = [];
  var outgoing_edges = [];

  var index = _.findIndex(nodes_array, function(o) { return o.id == params.nodes[0]; });
  node_label = nodes_array[index].label.replace("\n","")

  nodes_label = _.keyBy(nodes_array, 'id') 
  neighbor_nodes = [];

  _.forEach(network.getConnectedNodes(params.nodes[0]),function(i) {
    updateNode({id: i, color: {background: "#ffdb99"}})
    neighbor_nodes.push(i);
  })

  _.forEach(edges_array, function(i){
    if (i.from === params.nodes[0]) {
      outgoing_edges.push({id: i.to, label: i.label.replace("\n",""), node: nodes_label[i.to].label})
    }
    else if (i.to === params.nodes[0]) {
      incoming_edges.push({id: i.from, label: i.label.replace("\n",""), node: nodes_label[i.from].label})
    }
  })

  // node info
  out = ""
  out += "<div class='select_section'><div class='title'>Selected (Node/Type):</div>";
  out += "<div class='node_label'>" + node_label + "</div>"
  out += "</div>"

  // incoming edges
  if (incoming_edges.length > 0) {
    out += "<div class='select_section'>"
    out += "<div class='title'>Incoming edges [?s, ?p, :selection] (<span class='highlight'>" + incoming_edges.length + "</span>)</div>"
    out += "<table class='scroll'>"
    out += "<thead><tr><th>Subject Type</th><th>Predicate</th></tr></thead>"
    out += "<tbody>"
    _.forEach(incoming_edges,function(i){
      out += "<tr>"
      out += "<td><span style='color:blue; cursor:pointer' onclick='selectNode(["+i.id+"])'>" + i.node + "</a></td>"
      out += "<td>" + i.label + "</td>"
      out += "</tr>"
    })
    // out += "</ul>"
    out += "</tbody></table>"
    out += "</div>"
  }

  if (outgoing_edges.length > 0) {
    out += "<div class='select_section'>"
        out += "<div class='title'>Outgoing edges [:selection, ?p, ?o] (<span class='highlight'>" + outgoing_edges.length + "</span>)</div>"
    out += "<table class='scroll'>"
    out += "<thead><tr><th>Predicate</th><th>Object Type</th></tr></thead>"    
    out += "<tbody>"
    _.forEach(outgoing_edges,function(i){
      out += "<tr>"
      out += "<td>" + i.label + "</td>"
      out += "<td><span style='color:blue; cursor:pointer' onclick='selectNode(["+i.id+"])'>" + i.node + "</a></td>"
      out += "</tr>"
    })
    out += "</tbody></table>"
    out += "</div>"
  }

  return out;
}

function removeSelectedEdge(id) {
  removeEdge(id);
  $('#edge').html("");
}
function edgeInfo(params) {
  out = "<button onclick='removeSelectedEdge("+params.edges[0]+")'>Delete edge</button>"
  return out;
}

function draw() {

  // menu info
  graph_info = "<div class='title'> Graph ("
  graph_info += "<span class='highlight'>" + nodes_array.length + "</span> Types | "
  graph_info += "<span class='highlight'>" + edges_array.length + "</span> Relations)"
  graph_info += "<span style='float: right'><button onclick='reset()'>Reset</button></span></div>"
  $('#graph').html(graph_info)

  // create an array with nodes
  nodes = new vis.DataSet();
  nodes.add(nodes_array);

  // create an array with edges
  edges = new vis.DataSet();
  edges.add(edges_array);

  // create a network
  var container = document.getElementById('network');
  var data = {
    nodes: nodes,
    edges: edges
  };

  var options = {
    layout: {
      improvedLayout: false
    },
    interaction: {
      // multiselect: true
      hideEdgesOnDrag: true,
      // hideNodesOnDrag: true,
      // navigationButtons: true
    },
    physics:{
      barnesHut:{
        gravitationalConstant: -60000,
        springConstant: 0.03
      },
      timestep: 0.35,
      stabilization: {
        enabled:true,
        iterations:200,
        updateInterval:25
      }
    },
    nodes: {
      color: {
        border: "#003099",
        background: "#b2caff",
        highlight: {
          border: "red",
          background: "#ff9999",
        }
      }
    },
    edges: {
      //arrows: 'to'
    }
  };

  network = new vis.Network(container, data, options);

  network.on("selectNode", function (params) {
    // network.edges
    if (selected_node != params.nodes[0]){
      selected_node = params.nodes[0];
      $('#selection').html(nodeInfo(params));
      document.getElementById("menu").style.resize = "both";
      network.physics.stabilized = true;
      network.physics.options.enabled = false;
    }
  });
  network.on("dragStart",function(params){
    if (params.nodes.length>0 && selected_node != params.nodes[0]){
      if (selected_node != null) {
        _.forEach(neighbor_nodes,function(i) {
          updateNode({id: i, color: {background: "#b2caff"}})
        })
        neighbor_nodes = []
      }
      selected_node = params.nodes[0];
      $('#selection').html(nodeInfo(params));
      document.getElementById("menu").style.resize = "both";
      network.physics.stabilized = true;
      network.physics.options.enabled = false;
    }
  })
  network.on("deselectNode", function (params) {
    selected_node = null;
    $('#selection').html("");
    document.getElementById("menu").style.height = "";
    _.forEach(neighbor_nodes,function(i) {
      updateNode({id: i, color: {background: "#b2caff"}})
    })
    neighbor_nodes = [];
  });

  network.on("startStabilizing", function() {
    Pace.start()
  });
  network.once("stabilizationIterationsDone", function() {
    Pace.stop()
  });
}

// Slider to show/hide edges give their weight
function hideEdges(cutoff) {  
  for (var i=0; i<edges_array.length; i++) {
    if (edges_array[i].weight < cutoff) {
      edges_array[i].hidden = true;
      updateEdge(edges_array[i]);
    }
  }
}
function showEdges(cutoff) {
  for (var i=0; i<edges_array.length; i++) {
    if (edges_array[i].weight >= cutoff) {
      edges_array[i].hidden = false;
      updateEdge(edges_array[i]);
    }
  }
}
// also hide nodes without edges
function hideNodes(cutoff) {  
  for (var i=0; i<nodes_array.length; i++) {
    if (nodes_array[i].weight < cutoff) {
      nodes_array[i].hidden = true;
      updateNode(nodes_array[i]);
    }
  }
}
function showNodes(cutoff) {
  for (var i=0; i<nodes_array.length; i++) {
    if (nodes_array[i].weight >= cutoff) {
      nodes_array[i].hidden = false;
      updateNode(nodes_array[i]);
    }
  }
}

function percentile(arr, p) {
  if (arr.length === 0) return 0;
  if (typeof p !== 'number') throw new TypeError('p must be a number');
  if (p <= 0) return arr[0];
  if (p >= 1) return arr[arr.length - 1];

  var index = arr.length * p,
      lower = Math.floor(index),
      upper = lower + 1,
      weight = index % 1;

  if (upper >= arr.length) return arr[lower];
  return arr[lower] * (1 - weight) + arr[upper] * weight;
}

var edges_range = [];
edges_range_all = _.map(_.sortBy(edges_array,'weight'),'weight')
edges_range.push(percentile(edges_range_all,0))
edges_range.push(percentile(edges_range_all,0.1))
edges_range.push(percentile(edges_range_all,0.2))
edges_range.push(percentile(edges_range_all,0.3))
edges_range.push(percentile(edges_range_all,0.4))
edges_range.push(percentile(edges_range_all,0.5))
edges_range.push(percentile(edges_range_all,0.6))
edges_range.push(percentile(edges_range_all,0.7))
edges_range.push(percentile(edges_range_all,0.8))
edges_range.push(percentile(edges_range_all,0.9))
edges_range.push(percentile(edges_range_all,1))

const slider_node = document.querySelectorAll('input[id=node_slider]');
$(slider_node).rangeslider({
  polyfill: false,
  // Callback function
  onInit: function() {
    // console.log("init slider ")
  },
  // Callback function
  onSlide: function(position, value) {
    $('#cutoff_node').html(Math.round(edges_range[value]))
  },
  // // Callback function
  onSlideEnd: function(position, value) {
    network.physics.options.enabled = true;
    network.physics.stabilized = false;
    // need a split distribution of 10 cutoffs
    if (value < currentslider_node) {
      showEdges(edges_range[value])
      showNodes(edges_range[value])
    } else if (value > currentslider_node) {
      hideEdges(edges_range[value])
      hideNodes(edges_range[value])
    }
    currentslider_node = value;
  }
});

const slider_edge = document.querySelectorAll('input[id=edge_slider]');
$(slider_edge).rangeslider({
  polyfill: false,
  // Callback function
  onInit: function() {
    // console.log("init slider ")
  },
  // Callback function
  onSlide: function(position, value) {
    $('#cutoff_edge').html(Math.round(edges_range[value]))
  },
  // // Callback function
  onSlideEnd: function(position, value) {
    network.physics.options.enabled = true;
    // network.stabilize();
    Pace.start()
    // need a split distribution of 10 cutoffs
    if (value < currentslider_edge) {
      showEdges(edges_range[value])
      showNodes(edges_range[value])
    } else if (value > currentslider_edge) {
      hideEdges(edges_range[value])
      hideNodes(edges_range[value])
    }
    currentslider_edge = value;
    Pace.stop()
  }
});

function reset() {
  // location.reload();
  network.physics.options.enabled = true;
  network.physics.stabilize();
  network.fit()
}

// On Load
var currentslider_edge = 5;
var currentslider_node = 5;
draw()
hideEdges(edges_range[currentslider_edge])
hideNodes(edges_range[currentslider_edge])
$('#cutoff_edge').html(Math.round(edges_range[currentslider_edge]))
$('#cutoff_node').html(Math.round(edges_range[currentslider_edge]))
document.getElementById("menu").style.resize = "both";

function selectNode(id) {
  if (selected_node != null) {
    _.forEach(neighbor_nodes,function(i) {
      updateNode({id: i, color: {background: "#b2caff"}})
    })
    neighbor_nodes = [];
  }
  network.selectNodes([id]);
  params_tmp = { nodes: id }
  $('#selection').html(nodeInfo(params_tmp));
  document.getElementById("menu").style.resize = "both";
  network.physics.stabilized = true;
  network.physics.options.enabled = false;
  selected_node = id;
}
