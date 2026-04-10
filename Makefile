NAME = 109titration
CARGO =	cargo
PROFILE	= release
SRC_BIN	= target/$(PROFILE)/$(NAME)

.PHONY: all clean fclean re

all: $(NAME)

$(NAME):
	@$(CARGO) build --$(PROFILE) -q
	@cp $(SRC_BIN) ./$(NAME)

clean:
	@$(CARGO) clean -q

fclean: clean
	@rm -f ./$(NAME)

re: fclean all
