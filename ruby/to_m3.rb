# convert memory-1 or memory-2 strategy into memory-3 strategy

def m1_to_m3(m1)
  raise "invalid input" unless m1.size == 4 and m1.each_char.all? {|c| c=='c' or c=='d'}
  64.times.map do |i|
    b_1 = i&1
    a_1 = (i&8)/8
    idx = a_1*2 + b_1
    m1[idx]
  end.join
end

def m2_to_m3(m2)
  raise "invalid input" unless m1.size == 16 and m1.each_char.all? {|c| c=='c' or c=='d'}
  64.times.map do |i|
    b_1 = i&1
    b_2 = i&2/2
    a_1 = (i&8)/8
    a_2 = (i&16)/16
    idx = a_2*8 + a_1*4 + b_2*2 + b_1
    m1[idx]
  end.join
end

ARGV.each do |arg|
  if arg.size == 4
    puts m1_to_m3(arg)
  elsif arg.size == 16
    puts m2_to_m3(arg)
  else
    raise StandardError("unsupported format: #{arg}")
  end
end

