#https://stackoverflow.com/questions/9876227/rabbitmq-consume-one-message-if-exists-and-quit
# https://www.rabbitmq.com/confirms.html
# https://stackoverflow.com/questions/53426824/how-to-acknowledge-message-in-rabbitmq-in-a-different-channel
def receive():
    parameters = pika.ConnectionParameters(RabbitMQ_server)
    connection = pika.BlockingConnection(parameters)
    channel = connection.channel()
    channel.queue_declare(queue='toM')
    method_frame, header_frame, body = channel.basic_get(queue = 'toM')
    if method_frame.NAME == 'Basic.GetEmpty':
        connection.close()
        return ''
    else:
        channel.basic_ack(delivery_tag=method_frame.delivery_tag)
        connection.close()
        return body