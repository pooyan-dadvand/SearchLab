/* -*- c++ -*- */
#pragma once

class SLexception : public std::exception {
public:
  SLexception( const std::string msg) : m_Msg( msg) {}
  const char *what() const noexcept { return m_Msg.c_str(); }
private:
    std::string m_Msg;
};
