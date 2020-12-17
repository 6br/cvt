require "spec_helper"

RSpec.describe Conv do
  it "has a version number" do
    expect(Conv::VERSION).not_to be nil
  end

  it "should rehash without errors" do
    expect {Conv::CLI.new.invoke(:rehash, [], {})}.not_to raise_exception
  end

  it "should display output when no candidates" do
    expect {Conv::CLI.new.invoke(:convert, ["spec/sample/a.sam", "b.non"], {hidelist: true})}.to output("No candidates.").to_stdout
  end

  it "should display error when no file" do
    expect {Conv::CLI.new.invoke(:convert, ["spec/sample/no.sam", "b.bam"], {hidelist: true})}.to raise_exception(RuntimeError)
  end

  it "should display list" do
    expect {Conv::CLI.new.invoke(:list, ["spec/sample/a.sam"], {showonly: true, overwrite: true})}.to output("[\"bam\"] samtools view -bS spec/sample/a.sam > spec/sample/a.bam  (samtools)[\"gtf\"] cufflinks spec/sample/a.sam  (cufflinks)\n").to_stdout
  end

end

