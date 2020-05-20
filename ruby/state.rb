require 'pp'

class State

  A_STATES = [
      [:c,:c,:c],
      [:c,:c,:d],
      [:c,:d,:c],
      [:c,:d,:d],
      [:d,:c,:c],
      [:d,:c,:d],
      [:d,:d,:c],
      [:d,:d,:d]
  ].map {|x| x.freeze}.freeze

  B_STATES = Marshal.load( Marshal.dump(A_STATES) ).freeze

  ALL_STATES = A_STATES.product(B_STATES).map {|a,b| (a+b).freeze }.freeze

  def self.valid?(state)
    ALL_STATES.include?( state )
  end

  def self.index( state )
    ALL_STATES.index( state )
  end

  def self.make_from_id( id )
    raise "invalid arg: #{id}" if id < 0 or id > 63
    a = (0..5).reverse_each.map do |i|
      id[i] == 1 ? :d : :c
    end
    self.new(*a)
  end

  def self.make_from_str(s)
    raise unless s =~ /\A[cd]{6}\z/
    self.new( *s.each_char.map(&:to_sym) )
  end

  def self.c(str)
    str
      .gsub('c',"\e[42m\e[30mc\e[0m")
      .gsub('d',"\e[45m\e[30md\e[0m")
  end

  attr_reader :a_3,:a_2,:a_1,:b_3,:b_2,:b_1

  def initialize(a_3,a_2,a_1,b_3,b_2,b_1)
    @a_3 = a_3
    @a_2 = a_2
    @a_1 = a_1
    @b_3 = b_3
    @b_2 = b_2
    @b_1 = b_1
    unless [@a_3,@a_2,@a_1,@b_3,@b_2,@b_1].all? {|a| a == :d or a == :c }
      raise "invalid state"
    end
  end

  def to_a
    [@a_3,@a_2,@a_1,@b_3,@b_2,@b_1]
  end

  def to_s
    to_a.join('')
  end

  def inspect
    self.class.c(to_s)
  end

  def ==(other)
    self.to_id == other.to_id
  end

  def to_id
    to_a.map {|x| x==:d ? '1' : '0' }.join.to_i(2)
  end
  alias :to_i :to_id

  def next_state(act_a,act_b)
    self.class.new(@a_2,@a_1,act_a,@b_2,@b_1,act_b)
  end

  def prev_state(act_a,act_b)
    self.class.new(act_a,@a_3,@a_2,act_b,@b_3,@b_2)
  end

  def possible_prev_states
    [:c,:d].map do |a|
      [:c,:d].map do |b|
        prev_state(a,b)
      end
    end.flatten
  end

  def swap
    self.class.new(@b_3,@b_2,@b_1,@a_3,@a_2,@a_1)
  end

  def relative_payoff
    a = @a_1
    b = @b_1

    if a == b
      return 0
    elsif a == :c and b == :d
      return -1
    elsif a == :d and b == :c
      return 1
    else
      raise "must not happen"
    end
  end
end

if __FILE__ == $0
  require 'minitest/autorun'

  class StateTest < Minitest::Test

    def test_alld
      fs = State.make_from_id(63)
      assert_equal [:d,:d,:d,:d,:d,:d], fs.to_a
      assert_equal 63, fs.to_id
      assert_equal 0, fs.relative_payoff
    end

    def test_allc
      fs = State.make_from_id(0)
      assert_equal [:c,:c,:c,:c,:c,:c], fs.to_a
      assert_equal 0, fs.to_id
      assert_equal 0, fs.relative_payoff
    end

    def test_state43
      fs = State.make_from_id(43)
      assert_equal [:d, :c, :d, :c, :d, :d], fs.to_a
      assert_equal 43, fs.to_id
      assert_equal 0, fs.relative_payoff
      assert_equal [:c,:d,:c,:d,:d,:d], fs.next_state(:c,:d).to_a
      assert_equal 'cdcdcd', fs.prev_state(:c,:d).to_s
    end

    def test_state44
      fs = State.make_from_id(44)
      assert_equal [:d, :c, :d, :d, :c, :c], fs.to_a
      assert_equal 1, fs.relative_payoff
      ns = fs.next_state(:c,:d)
      assert_equal [:c,:d,:c,:c,:c,:d], ns.to_a
      assert_equal -1, ns.relative_payoff
      assert_equal 'ddccdc', fs.prev_state(:d,:c).to_s
    end

    def test_swap
      s = State.make_from_id(46)
      assert_equal [:d, :c, :d, :d, :d, :c], s.to_a
      s2 = s.swap
      assert_equal [:d, :d, :c, :d, :c, :d], s2.to_a
    end

    def test_equality
      fs1 = State.make_from_id(15)
      fs2 = State.new(:c,:c,:d,:d,:d,:d)
      assert_equal true, fs1 == fs2
    end

    def test_make_from_str
      s = State.make_from_str('dcdccd')
      assert_equal State.new(:d,:c,:d,:c,:c,:d), s
    end

    def test_prev_states
      s = State.make_from_str('dcdccd')
      a = s.possible_prev_states.map(&:to_s)
      assert_equal ['cdcccc','cdcdcc','ddcccc','ddcdcc'], a
    end
  end
end

