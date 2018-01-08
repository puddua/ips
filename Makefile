DOC= doxygen
SUBDIRS = src test

subdirs:
	for dir in $(SUBDIRS); do \
	$(MAKE) -C $$dir; \
	done

doc:
	$(DOC) 

clean:
	for dir in $(SUBDIRS); do \
	$(MAKE) clean -C $$dir; \
	done
	rm -r html/
