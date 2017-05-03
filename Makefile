
all : run-without run-with

run-without:
	@python without.py

run-with:
	@python with.py

distclean : clean uninstall

clean :
	@rm -f *~

uninstall:
	@rm -rf __pycache__


