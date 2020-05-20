require_relative 'state'
require_relative 'strategy'
require_relative 'graph'

module DFAminimize

  class UnionFind
    def initialize(n)
      @par  = Array.new(n) {|i| i}
    end

    def root(i)
      if @par[i] != i
        r = root(@par[i])
        @par[i] = r
      end
      @par[i]
    end

    def roots
      to_h.keys.sort
    end

    def merge(i, j)
      ri = root(i)
      rj = root(j)
      return false if ri == rj  # already merged
      ri,rj = rj,ri if ri > rj
      @par[rj] = ri
    end

    def to_h
      @par.size.times.group_by {|i| root(i) }
    end
  end

  def self._equivalent(str, i, j, uf, noisy)
    # both action & next state must be identical
    raise unless str.action(i) == str.action(j)
    act_a = str.action(i)
    err_a = act_a == :c ? :d : :c
    acts_other = [:c,:d]
    acts_other.each do |act_b|
      ni = State.make_from_id(i).next_state(act_a, act_b).to_id
      nj = State.make_from_id(j).next_state(act_a, act_b).to_id
      return false unless uf.root(ni) == uf.root(nj)
      if noisy
        ni2 = State.make_from_id(i).next_state(err_a, act_b).to_id
        nj2 = State.make_from_id(j).next_state(err_a, act_b).to_id
        return false unless uf.root(ni2) == uf.root(nj2)
      end
    end
    return true
  end

  def self.minimize_DFA(str, noisy=false)
    uf_0 = UnionFind.new(Strategy::N)
    # initial grouping by the action c/d
    c_rep = Strategy::N.times.find {|i| str.action(i) == :c}
    d_rep = Strategy::N.times.find {|i| str.action(i) == :d}
    Strategy::N.times do |i|
      t = (str.action(i)==:c ? c_rep : d_rep)
      uf_0.merge(i, t)
    end

    loop do
      uf_0_h = uf_0.to_h
      uf = UnionFind.new(Strategy::N)
      uf_0_h.each do |r,s|  # refinint a set in uf_0
        s.combination(2).each do |i,j|
          if _equivalent(str, i, j, uf_0, noisy)
            uf.merge(i,j)
          end
        end
      end
      break if uf.to_h == uf_0_h
      uf_0 = uf
    end

    g = DirectedGraph.new(Strategy::N)
    org_g = noisy ? str.noisy_transition_graph : str.transition_graph
    org_g.for_each_link do |i,j|
      ri = uf_0.root(i)
      rj = uf_0.root(j)
      g.add_link(ri,rj)
    end
    g.remove_duplicate_links!

    return uf_0, g
  end
end

if __FILE__ == $0 and ARGV.size == 1
  DEBUG = true

  s = ARGV[0]
  raise "unsupported input format" unless s.length == 64

  s = s.gsub('*', 'c')
  str = Strategy.make_from_str(s)
  $stderr.puts str.inspect
  $stderr.puts "defensible?      : #{str.defensible?}"
  $stderr.puts "efficient?       : #{str.efficient?}"
  $stderr.puts "distinguishable? : #{str.distinguishable?}"

  # automaton representation
  full_rep = true
  uf, min_g = DFAminimize.minimize_DFA(str, full_rep)
  $stderr.puts uf.to_h.inspect
  $stderr.puts "automaton size : #{uf.to_h.size}"

  def trace_path(str, init_state = 'cccccd')
    path = []
    s = State.make_from_str(init_state)
    until path.include?(s)
      path.push(s)
      s = str.next_state_with_self(s)
    end
    path
  end

  path = trace_path(str)
  $stderr.puts "path from 'ccc,ccd'"
  $stderr.puts "  recovered in #{path.length-1} rounds"
  $stderr.puts path.map {|s| "#{s} (#{uf.root(s.to_id)},#{uf.root(s.swap.to_id)})" }.join(' -> ')

  path = trace_path(str, 'dddddc')
  $stderr.puts "path from 'ddd,ddc'"
  $stderr.puts "  ends in #{path.length-1} rounds"
  $stderr.puts path.map {|s| "#{s} (#{uf.root(s.to_id)},#{uf.root(s.swap.to_id)})" }.join(' -> ')

  def to_dot(str, uf, min_g, full_rep)
    mapped = uf.roots.map do |n|
      c = (str.action(n) == :c) ? "#afeeee" : "#ffa07a"
      [ n, {label: "#{str.action(n)}@#{n}", style: "filled", fillcolor: c} ]
    end
    attr = Hash[mapped]
    link_label = Hash.new {|h,k| h[k] = [] }
    t = uf.to_h
    g = full_rep ? str.noisy_transition_graph : str.transition_graph
    g.for_each_link do |i,j|
      e = [uf.root(i), uf.root(j)]
      last_ac = [8,1].map{|m| ((j&m)==m)?'d':'c'}.join
      last_ac ='*'+last_ac if full_rep and last_ac[0] != str.action(i).to_s
      link_label[e].push( last_ac )
    end
    link_label = link_label.map {|k,v| [k, v.sort.uniq.join(',')] }.to_h
    sio = StringIO.new
    min_g.to_dot(remove_isolated: true, node_attributes: attr, edge_labels: link_label)
  end

  File.open('a.dot', 'w') do |io|
    io.puts to_dot(str, uf, min_g, full_rep)
  end
  $stderr.puts "a.dot was written"
end

