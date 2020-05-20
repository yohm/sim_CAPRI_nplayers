require 'pp'

class DirectedGraph

  attr_reader :n, :links

  def initialize(size)
    @n = size
    @links = Hash.new {|hsh,key| hsh[key] = Array.new }
  end

  def add_link( from, to )
    raise "invalid index: #{from}->#{to}" unless from < @n and to < @n
    @links[from].push(to)
  end

  def sccs
    f = ComponentFinder.new(self)
    f.strongly_connected_components
  end

  def terminanl_components
    s = sccs
    has_self_loops = s.select {|scc| scc.size == 1}.flatten.select do |n|
      @links[n].include?(n)
    end
    s.select {|scc| scc.size > 1} + has_self_loops.map {|n| [n]}
  end

  def transient_nodes
    isolated_nodes = sccs.select {|scc| scc.size == 1 }.flatten
    has_selfloop = isolated_nodes.select do |n|
      @links[n].include?(n)
    end
    isolated_nodes - has_selfloop
  end

  def non_transient_nodes
    (0..(@n-1)).to_a - transient_nodes
  end

  def remove_duplicate_links!
    @links.values.each {|ns| ns.uniq! }
  end

  # node_attributes = { 0=> {label: "0_cccccc", fontcolor: "red"}, 1=>{ ...}, ... }
  # edge_attributes = { [0,1] => {color: "yellow"}, [1,3]=>{....}, ...}
  # for available attributes, see the graphviz documentation https://graphviz.gitlab.io/_pages/doc/info/attrs.html#d:fillcolor
  def to_dot(node_attributes: {}, remove_isolated: false, node_ranks: [], edge_labels: {})
    io = StringIO.new
    io.puts "digraph \"\" {"
    @n.times do |ni|
      next if remove_isolated and @links[ni].empty?
      label = node_attributes.dig(ni,:label) || ni.to_s
      fontcolor = node_attributes.dig(ni,:fontcolor) || "black"
      a = {label: label, fontcolor: fontcolor}.merge( node_attributes.dig(ni) || {} )
      s = a.each_pair.map {|k,v| "#{k}=\"#{v}\"" }.join('; ')
      io.puts "  #{ni} [ #{s} ];"
      # io.puts "  #{ni} [ label=\"#{label}\"; fontcolor = #{fontcolor} ];"
    end
    @n.times do |ni|
      next if remove_isolated and @links[ni].empty?
      @links[ni].each do |nj|
        io.puts "  #{ni} -> #{nj} [label=\"#{edge_labels[[ni,nj]]}\"];"
      end
    end
    if node_ranks.size > 0
      ranks = node_ranks.map.with_index {|_,i| "rank_#{i}"}
      io.puts "  #{ranks.join(' -> ')}"
      node_ranks.each_with_index do |nodes,i|
        io.puts "  {rank=same; #{ranks[i]}; #{nodes.join(';')};}"
      end
    end
    io.puts "}"
    io.flush
    io.string
  end

  def is_accessible?(from, to)
    found = false
    dfs(from) {|n|
      found = true if n == to
    }
    found
  end

  def for_each_link
    @n.times do |ni|
      @links[ni].each do |nj|
        yield ni, nj
      end
    end
  end

  def dfs(start, &block)
    stack=[]
    dfs_impl = lambda do |n|
      block.call(n)
      stack.push(n)
      @links[n].each do |nj|
        next if stack.include?(nj)
        dfs_impl.call(nj)
      end
    end
    dfs_impl.call(start)
  end

  def self.common_subgraph(g1,g2)
    g = self.new( g1.n )
    links1 = []
    g1.for_each_link {|ni,nj| links1.push( [ni,nj] ) }
    links2 = []
    g2.for_each_link {|ni,nj| links2.push( [ni,nj] ) }
    common_links = links1 & links2
    common_links.each {|l| g.add_link(*l) }
    g
  end
end

class DirectedWeightedGraph < DirectedGraph

  attr_reader :weights

  def initialize(size)
    super size
    @weights = Hash.new {|h,k| h[k] = {} }
  end

  def add_link( from, to, weight )
    super(from, to)
    @weights[from][to] = weight
  end

  def bellman_ford(from)
    dist = Array.new(@n) {|i| i==from ? 0 : Float::INFINITY }
    (@n-1).times do |t|
      for_each_link do |u,v,w|
        dist[v] = dist[u] + w if dist[v] > dist[u] + w
      end
    end
    # check negative cycle
    for_each_link do |u,v,w|
      raise "negative cycle detected" if dist[u] + w < dist[v]
    end
    dist
  end

  def for_each_link
    super do |i,j|
      yield i, j, @weights[i][j]
    end
  end

  def has_negative_cycle?
    d = Array.new(@n) {|i| Array.new(@n,Float::INFINITY) }
    for_each_link do |i,j,w|
      d[i][j] = w
    end
    @n.times do |i|
      d[i][i] = 0 if d[i][i] > 0
    end
    @n.times do |k|
      @n.times do |i|
        @n.times do |j|
          x = d[i][k] + d[k][j]
          d[i][j] = x if d[i][j] > x
          if j == i and d[i][j] < 0
            return true
          end
        end
      end
    end
    false
  end
end

class ComponentFinder

  def initialize( graph )
    @g = graph
    @t = 0

    @desc =  Array.new(@g.n, nil)
    @low  =  Array.new(@g.n, nil)
    @stack = []
    @on_stack  =  Array.new(@g.n, false)
  end

  def strongly_connected_components
    @sccs = []
    @g.n.times do |v|
      if @desc[v].nil?
        strong_connect(v)
      end
    end
    @sccs
  end

  private
  def strong_connect(v)
    @desc[v] = @t
    @low[v] = @t
    @t += 1

    @stack.push(v)
    @on_stack[v] = true

    @g.links[v].each do |w|
      if @desc[w].nil?
        strong_connect(w)
        @low[v] = @low[w] if @low[w] < @low[v]
      elsif @on_stack[w]
        @low[v] = @desc[w] if @desc[w] < @low[v]
      end
    end

    # if v is a root node, pop the stack and generate an scc
    scc = []
    if @low[v] == @desc[v]
      loop do
        w = @stack.pop
        @on_stack[w] = false
        scc.push(w)
        break if v == w
      end
      @sccs.push( scc )
    end
  end
end

if __FILE__ == $0
  require 'minitest/autorun'

  class TestDirectedGraph < Minitest::Test

    def setup
      @g = DirectedGraph.new(5)
      @g.add_link(1, 0)
      @g.add_link(0, 2)
      @g.add_link(2, 1)
      @g.add_link(0, 3)
      @g.add_link(3, 4)
      @g.add_link(4, 4)
    end

    def test_add_link
      assert_equal 5, @g.n
      expected = {0=>[2,3],1=>[0],2=>[1],3=>[4],4=>[4]}
      assert_equal expected, @g.links
    end

    def test_sccs
      assert_equal [ [0,1,2], [3], [4] ], @g.sccs.map(&:sort).sort
    end

    def test_transient_nodes
      assert_equal [3], @g.transient_nodes
      assert_equal [0,1,2,4], @g.non_transient_nodes.sort
      assert_equal @g.terminanl_components.flatten.sort, @g.non_transient_nodes.sort
    end

    def test_to_dot
      sio = StringIO.new
      dot = @g.to_dot
      assert_equal dot.empty?, false
    end

    def test_dfs
      traversed = []
      @g.dfs(0) {|n| traversed << n }
      assert_equal [0,2,1,3,4], traversed
    end

    def test_accessible
      assert_equal true, @g.is_accessible?(0,4)
      assert_equal false, @g.is_accessible?(3,0)
    end
  end

  class TestDirectedWeightedGraph < Minitest::Test

    def setup
      @g = DirectedWeightedGraph.new(5)
      @g.add_link(1, 0, 0.1)
      @g.add_link(0, 2, 0.2)
      @g.add_link(2, 1, 0.3)
      @g.add_link(0, 3, 0.4)
      @g.add_link(3, 4, 0.5)
      @g.add_link(4, 4, 0.6)
    end

    def test_add_link
      assert_equal 5, @g.n
      expected = {0=>[2,3],1=>[0],2=>[1],3=>[4],4=>[4]}
      assert_equal expected, @g.links
    end

    def test_weights
      w_expected = {0=>{2=>0.2,3=>0.4} ,1=>{0=>0.1}, 2=>{1=>0.3}, 3=>{4=>0.5}, 4=>{4=>0.6} }
      assert_equal w_expected, @g.weights
    end

    def test_for_each_link
      count = 0
      @g.for_each_link do |i,j,w|
        count += 1
        assert_equal @g.weights[i][j], w
      end
      assert_equal count, 6
    end

    def test_bellman_ford
      dist = @g.bellman_ford(0)
      expected = [0, 0.5, 0.2, 0.4, 0.9]
      assert_equal dist, expected
    end

    def test_negative_cycle
      assert_equal false, @g.has_negative_cycle?
      @g.add_link( 0, 0, -1 )
      assert_equal true, @g.has_negative_cycle?
    end

    def test_negative_cycle2
      g = DirectedWeightedGraph.new(3)
      g.add_link(0, 1, 1)
      g.add_link(1, 2, 1)
      g.add_link(2, 0, -3)
      assert_equal true, g.has_negative_cycle?
    end
  end
end

